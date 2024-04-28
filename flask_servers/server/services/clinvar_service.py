import requests
import xmltodict
from flask import current_app
import pandas as pd
import time
import json
from Bio import Entrez
from requests import RequestException

from typing import Dict, List
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
def retrieve_clinvar_dict(clinvar_uid: str):
    # retrieve clinvar info for a single variant
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id={clinvar_uid}&from_esearch=true"

    try:
        clinvar_res = requests.post(url)
    except RequestException as e:
        current_app.logger.error(f'Failed to connect to Entrez Clinvar Service: {e}')
        # send service unavailable status code
        return InternalResponse(None, 503, e)

    if clinvar_res.status_code != 200:
        current_app.logger.error(f'Entrez Clinvar Service failed: {clinvar_res.reason}')
        return InternalResponse(None, clinvar_res.status_code, clinvar_res.reason)
    else:
        clinvar_res_dict = xmltodict.parse(clinvar_res.text)
        return InternalResponse(clinvar_res_dict, 200)


def retrieve_multiple_clinvar_dict(clinvar_ids: List[str]): # TODO: use this to replace retrieve_clinvar_dict
    # retrieve clinvar info for multiple variants
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id={','.join(clinvar_ids)}&from_esearch=true"

    try:
        clinvar_res = requests.post(url)
    except RequestException as e:
        current_app.logger.error(f'Failed to connect to Entrez Clinvar Service: {e}')
        # send service unavailable status code
        return InternalResponse(None, 503, e)

    if clinvar_res.status_code != 200:
        current_app.logger.error(f'Entrez Clinvar Service failed: {clinvar_res.reason}')
        return InternalResponse(None, clinvar_res.status_code, clinvar_res.reason)
    else:
        clinvar_res_dict = xmltodict.parse(clinvar_res.text)
        return InternalResponse(clinvar_res_dict, 200)


# Compare variant's properties with the properties of the ClinVar variant.
# Note:
# - When a ClinVar variant has multiple genes, then each gene is compared with the expeceted gene until a match is found.
# - The genotype is compared when the Variant Validator request is made (not in the below fn).
def compare_clinvar_variant_with_expected_variant(genome_version: str, retrieved_var_clinvar_dict, gene: str,
                                                  chr: str, chr_pos: str) -> tuple[bool, str]:

    clinvar_allele = (retrieved_var_clinvar_dict.get('ClinVarResult-Set').get('VariationArchive')
                      .get('ClassifiedRecord').get('SimpleAllele'))

    # compare expected and retrieved variant's gene
    clinvar_genes = clinvar_allele.get('GeneList')

    # if none of the genes match the expected gene
    if gene not in [clinvar_genes[key]['@Symbol'] for key in clinvar_genes.keys()]:
        return False, (f"None of the gene names {[clinvar_genes[key]['@Symbol'] for key in clinvar_genes.keys()]} "
                       f"match the expected gene name {gene}!")

    clinvar_locations = clinvar_allele.get('Location').get('SequenceLocation')
    for loc in clinvar_locations:
        if loc.get('@Assembly') == genome_version:
            # compare expected and retrieved variant's chromosome
            if loc.get('@Chr') != chr:
                return False, f"Chromosome {loc.get('@Chr')} does not match the expected chromosome {chr}!"
            # compare expected and retrieved variant's chromosome position
            elif loc.get('@start') != str(chr_pos):
                return False, (f"Chromosome start position {loc.get('@start')} does not match the expected chromosome "
                               f"position {chr_pos}!")
            break
        #TODO compare ref & alt alleles

    return True, ''


# Retrieve the ClinVar variant's clinical significance.
def extract_clinvar_germline_classification(clinvar_dict: Dict):
    germline_classification = (clinvar_dict.get('ClinVarResult-Set').get('VariationArchive').get('ClassifiedRecord')
                      .get('Classifications').get('GermlineClassification'))

    last_eval = ""
    if germline_classification.get('@DateLastEvaluated') != "1/01/01 00:00":
        last_eval = germline_classification.get('@DateLastEvaluated').replace('-', '/') + ' 00:00'

    return {'description': germline_classification.get('Description'), 'last_evaluated': last_eval,
            'review_status': germline_classification.get('ReviewStatus')}


# Retrieve the ClinVar variant's canonical SPDI.
def extract_clinvar_canonical_spdi(clinvar_dict: Dict):
    canonical_spdi = (clinvar_dict.get('ClinVarResult-Set').get('VariationArchive').get('ClassifiedRecord')
                         .get('SimpleAllele').get('CanonicalSPDI'))

    return canonical_spdi


# Retrieve the ClinVar variant's uid.
def extract_clinvar_uid(clinvar_dict: Dict):
    return clinvar_dict.get('ClinVarResult-Set').get('VariationArchive').get('@VariationID')


def clinvar_clinical_significance_pipeline(genome_version: str, rsid: str, gene: str, chr: str,
                                           chr_pos: str) -> InternalResponse:
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
        var_clinvar_ids = retrieve_clinvar_ids_res.data

        # check the number of ClinVar IDs returned for a variant search
        if len(var_clinvar_ids) == 0:
            error_msg = 'ClinVar ID has not been found!'
            is_success = False
        elif len(var_clinvar_ids) > 1:
            error_msg = f'The following ClinVar IDs returned: {var_clinvar_ids}'
            is_success = False
        # if only a single ClinVar ID has been returned
        else:
            var_clinvar_id = var_clinvar_ids[0]

            time.sleep(0.5)

            retrieve_clinvar_dict_res = retrieve_clinvar_dict(var_clinvar_id)

            if retrieve_clinvar_dict_res.status != 200:
                current_app.logger.error(
                    f"Retrieval of ClinVar document summary for ClinVar Id {var_clinvar_id} and RSID {rsid} failed 500!")
                return InternalResponse(None, 500)
            else:
                clinvar_dict = retrieve_clinvar_dict_res.data

                if is_success:
                    are_equivalent, error_msg = compare_clinvar_variant_with_expected_variant(genome_version,
                                                                                              clinvar_dict,
                                                                                              gene, chr, chr_pos)
                    if are_equivalent:
                        clinical_significance = extract_clinvar_germline_classification(clinvar_dict)
                        canonical_spdi = extract_clinvar_canonical_spdi(clinvar_dict)
                        uid = extract_clinvar_uid(clinvar_dict)
                    else:
                        is_success = False

        return InternalResponse((is_success, clinical_significance, canonical_spdi, uid, error_msg), 200)


# TODO: retrieve other vital infor from clinvar (such as last modfied)
# Retrieve Clinvar variant classifications for every variant attempt to retrieve a corresponding
# ClinVar variant and extract its clinical significance.
def retrieve_clinvar_variant_classifications(vus_df: pd.DataFrame) -> InternalResponse:
    genome_version = 'GRCh37'
    # performance_dict = {}

    # make a copy of the dataframe to be able to iterate through it whilst modifying the original dataframe
    new_vus_df = vus_df.copy()

    for index, row in new_vus_df.iterrows():
        current_app.logger.info(
            f"Retrieving information for:\n\tGene: {row['Gene']}\n\tChromosome: {row['Chr']}\n\tChromosome position: "
            f"{row['Position']}\n\tGenotype: {row['Genotype']}")

        # TODO: add fix for when multiple rsids are found
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
