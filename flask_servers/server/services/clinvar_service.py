from flask import current_app, Response
import pandas as pd
import time
import json
from Bio import Entrez

from server.responses.internal_response import InternalResponse


# Function to retrieve the variant's ClinVar ID The Bio Entrez package is used to extract information from ClinVar.
# Entrez is name of the NCBI infrastructure which provides access to all of the NCBI (US) databases. This package is
# used to search for variants using the RSID as the search term. The esearch() function is used to retrieve unique
# identifiers from ClinVar for these variants. It is assumed that the first returned identifier is the most relevant
# one.
def retrieve_clinvar_ids(rsid: str) -> InternalResponse:
    try:
        # retrieve unique identifiers of variants
        var_ids_handle = Entrez.esearch(db="clinvar", term=rsid, retmax=1)
        var_ids_record = Entrez.read(var_ids_handle)
        var_ids_handle.close()

        ids = var_ids_record['IdList']
    except IOError as e:
        current_app.logger.error(f'Network error when calling Entrez.esearch(): {str(e)}')
        return InternalResponse(None, e.errno, str(e))

    return InternalResponse(ids, 200)


# Function to retrieve the variant's Clinvar document summary
# The esummary() function is used to return document summaries from ClinVar for a given ClinVar unique variant id.
def retrieve_clinvar_document_summary(clinvar_id: str):
    try:
        # retrieve document summary - can take multiple ids
        handle = Entrez.esummary(db="clinvar", id=clinvar_id, retmode='json')

        # Read the JSON data
        json_data = handle.read()

        # Parse the JSON data
        document_summary_dict = json.loads(json_data)
        handle.close()
    except IOError as e:
        current_app.logger.error(f'Network error when calling Entrez.esummary(): {str(e)}')
        return InternalResponse(None, e.errno, str(e))

    return InternalResponse(document_summary_dict['result'][clinvar_id], 200)


# Compare variant's properties with the properties of the ClinVar variant.
# Note:
# - When a ClinVar variant has multiple genes, then each gene is compared with the expeceted gene until a match is found.
# - The genotype is compared when the Variant Validator request is made (not in the below fn).
def compare_clinvar_variant_with_expected_variant(genome_version: str, retrieved_var_clinvar_doc_summary, gene: str,
                                                  chr: str, chr_pos: str) -> tuple[bool, str]:
    # compare expected and retrieved variant's gene
    doc_genes = retrieved_var_clinvar_doc_summary.get('genes')

    # if none of the genes match the expected gene
    if gene not in [x['symbol'] for x in doc_genes]:
        return False, f"None of the gene names {[gene_name['symbol'] for gene_name in doc_genes]} match the expected gene name {gene}!"

    variation_set_arr = retrieved_var_clinvar_doc_summary.get('variation_set')
    if len(variation_set_arr) > 0:
        variation_set = variation_set_arr[0]
        variation_locations = variation_set.get('variation_loc')

        for loc in variation_locations:
            if loc['assembly_name'] == genome_version:
                # compare expected and retrieved variant's chromosome
                if loc['chr'] != chr:
                    return False, f"Chromosome {loc['chr']} does not match the expected chromosome {chr}!"
                # compare expected and retrieved variant's chromosome position
                elif loc['start'] != str(chr_pos):
                    return False, f"Chromosome start position {loc['start']} does not match the expected chromosome position {chr_pos}!"
                break

    return True, ''


# Retrieve the ClinVar variant's clinical significance.
def extract_clinvar_clinical_significance(clinvar_doc_summary):
    clinical_significance_obj = clinvar_doc_summary.get('clinical_impact_classification')

    last_eval = ""
    if (clinical_significance_obj['last_evaluated'] != "1/01/01 00:00" and
            len(clinical_significance_obj['description']) == 0):
        last_eval = clinical_significance_obj['last_evaluated']

    return {'description': clinical_significance_obj['description'], 'last_evaluated': last_eval,
            'review_status': clinical_significance_obj['review_status']}


# Retrieve the ClinVar variant's canonical SPDI.
def extract_clinvar_canonical_spdi(clinvar_doc_summary):
    canonical_spdi = ''

    variation_set_arr = clinvar_doc_summary.get('variation_set')

    if len(variation_set_arr) > 0:
        canonical_spdi = variation_set_arr[0]['canonical_spdi']

    return canonical_spdi


# Retrieve the ClinVar variant's uid.
def extract_clinvar_uid(clinvar_doc_summary):
    return clinvar_doc_summary['uid']


def clinvar_clinical_significance_pipeline(genome_version: str, rsid: str, gene: str, chr: str, chr_pos: str) -> InternalResponse:
    is_success = True
    clinical_significance = {}
    canonical_spdi = ''
    uid = ''
    error_msg = ''

    retrieve_clinvar_ids_res = retrieve_clinvar_ids(rsid)

    if retrieve_clinvar_ids_res.status != 200:
        current_app.logger.error(f"Retrieval of ClinVar id for {rsid} failed 500!")
        return InternalResponse(None, 500)
    else:
        var_ids = retrieve_clinvar_ids_res.data

        # check the number of ClinVar IDs returned for a variant search
        if len(var_ids) == 0:
            error_msg = 'ClinVar ID has not been found!'
            is_success = False
        elif len(var_ids) > 1:
            error_msg = f'The following ClinVar IDs returned: {var_ids}'
            is_success = False
        # if only a single ClinVar ID has been returned
        else:
            var_id = var_ids[0]

            time.sleep(0.5)

            retrieve_clinvar_document_summary_res = retrieve_clinvar_document_summary(var_id)

            if retrieve_clinvar_document_summary_res.status != 200:
                current_app.logger.error(f"Retrieval of ClinVar document summary for ClinVar Id {var_id} and RSID {rsid} failed 500!")
                return InternalResponse(None, 500)
            else:
                document_summary = retrieve_clinvar_document_summary_res.data

                if is_success:
                    are_equivalent, error_msg = compare_clinvar_variant_with_expected_variant(genome_version, document_summary,
                                                                                              gene, chr, chr_pos)
                    if are_equivalent:
                        clinical_significance = extract_clinvar_clinical_significance(document_summary)
                        canonical_spdi = extract_clinvar_canonical_spdi(document_summary)
                        uid = extract_clinvar_uid(document_summary)
                    else:
                        is_success = False

        return InternalResponse((is_success, clinical_significance, canonical_spdi, uid, error_msg), 200)


# TODO: retrieve other vital infor from clinvar (such as last modfied)
# Retrieve Clinvar variant classifications for every variant attempt to retrieve a corresponding ClinVar variant and extract its clinical significance.
def retrieve_clinvar_variant_classifications(vus_df: pd.DataFrame) -> InternalResponse:
    genome_version = 'GRCh37'
    performance_dict = {}

    # make a copy of the dataframe to be able to iterate through it whilst modifying the original dataframe
    new_vus_df = vus_df.copy()

    for index, row in new_vus_df.iterrows():
        current_app.logger.info(
            f"Retrieving information for:\n\tGene: {row['Gene']}\n\tChromosome: {row['Chr']}\n\tChromosome position: {row['Position']}\n\tGenotype: {row['Genotype']}")

        #TODO: add fix for when multiple rsids are found
        clinvar_clinical_significance_pipeline_res = clinvar_clinical_significance_pipeline(genome_version,
                                                                                             row['RSID'], row['Gene'],
                                                                                             row['Chr'],
                                                                                             row['Position'])

        if clinvar_clinical_significance_pipeline_res.status != 200:
            current_app.logger.error(
                f"ClinVar clinical significance pipeline failed for variant with RSID {row['RSID']}!")
            return InternalResponse(None, 500)
        else:
            # execute pipeline
            is_success, clinical_significance, canonical_spdi, uid, error_msg = (
                clinvar_clinical_significance_pipeline_res.data)

            if clinical_significance:
                vus_df.at[index, 'Clinvar classification'] = clinical_significance['description']
                vus_df.at[index, 'Clinvar classification last eval'] = clinical_significance['last_evaluated']
                vus_df.at[index, 'Clinvar classification review status'] = clinical_significance['review_status']

            vus_df.at[index, 'Clinvar error msg'] = error_msg
            vus_df.at[index, 'Clinvar canonical spdi'] = canonical_spdi
            vus_df.at[index, 'Clinvar uid'] = uid

            # update performance for given variant
    #         if row['VUS Id'] in performance_dict:
    #             if is_success:
    #                 performance_dict[row['VUS Id']] = True
    #         else:
    #             performance_dict[row['VUS Id']] = is_success
    #
    # num_of_variants = len(performance_dict.keys())
    # success = len([key for key in performance_dict.keys() if performance_dict[key] is True])
    # failure = num_of_variants - success
    #
    # current_app.logger.info(
    #     f'Found in ClinVar: {round(success / num_of_variants * 100, 2)}%\nNot found in ClinVar: {round(failure / num_of_variants * 100, 2)}%')

    return InternalResponse(vus_df, 200)
