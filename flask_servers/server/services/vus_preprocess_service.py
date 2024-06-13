from datetime import datetime
from typing import List, Hashable, Dict

from flask import current_app, Response
import pandas as pd
from Bio import Entrez
from flask_login import current_user
from sqlalchemy.exc import SQLAlchemyError
from werkzeug.datastructures import FileStorage
import json
import re

from server import db
from server.helpers.data_helper import prep_vus_df_for_react, convert_df_to_list, prep_unprocessed_vus_dict_for_react
from server.helpers.db_access_helper import get_variant_from_db
from server.models import Variants, GeneAnnotations, GeneAttributes, DbSnp, \
    Clinvar, ExternalReferences, FileUploads, Samples, VariantsSamples, Genotype, \
    VariantsSamplesUploads, ManualUploads, AcmgRules, VariantsAcmgRules, VariantHgvs, Classification, Reviews
from server.responses.internal_response import InternalResponse
from server.services.dbsnp_service import get_rsids_from_dbsnp

from server.services.clinvar_service import retrieve_clinvar_variant_classifications, get_updated_external_references_for_existing_vus, store_clinvar_info
from server.services.phenotype_service import get_hpo_term_from_phenotype_name, append_phenotype_to_sample
from server.services.view_vus_service import get_last_saved_clinvar_update
from server.services.publications_service import retrieve_and_store_variant_publications

Entrez.email = "esther.spiteri.18@um.edu.mt"


# TODO: fn can be eliminated later on since user will be asked which of the genes he/she would like to use

# generate separate records for each gene when a VUS has multiple genes
def handle_vus_with_multiple_genes(vus_df: pd.DataFrame):
    multiple_genes = []
    rows_to_delete = []
    new_rows_dictionaries = []
    var_id = 0

    # insert id column so that variants with multiple genes have multiple entries in the dataframe with the same id
    vus_df.insert(0, 'VUS Id', 0)

    # make a copy of the dataframe to be able to iterate through it whilst modifying the original dataframe
    new_vus_df = vus_df.copy()

    # note VUS with multiple genes and generate records for each gene (same VUS ID)
    for index, row in new_vus_df.iterrows():
        # set variant id
        vus_df.at[index, 'VUS Id'] = var_id

        # get the gene name
        gene_string = row['Gene']
        genes = gene_string.split(',')

        # note if variant has multiple genes
        filtered_genes = [g for g in genes if 'AS1' not in g and 'LOC' not in g]
        if len(filtered_genes) > 1:
            multiple_genes.append(var_id)

        if len(genes) > 1:
            rows_to_delete.append(index)

            # create a variant record for every gene
            for i, gene in enumerate(genes):
                # exclude anti-sense and locus
                if 'AS1' not in gene and 'LOC' not in gene:
                    new_row = row.copy().to_dict()
                    new_row['Gene'] = gene
                    new_row['VUS Id'] = var_id
                    new_rows_dictionaries.append(new_row)

        # increment variant id
        var_id += 1

        # TODO: handle cases where not SNVs

    current_app.logger.info(f"Percentage of VUS with multiple genes: {round(len(multiple_genes) / var_id * 100, 2)}%\n"
                            f"VUS IDs of VUS with multiple genes: {[g for g in multiple_genes]}")

    # clean up - remove multiple gene records (since they will be replaced by individual records)
    # drop rows containing variants with multiple genes (including those with AS1 and LOC)
    vus_df = vus_df.drop([g for g in rows_to_delete])
    vus_df = vus_df.reset_index(drop=True)

    # insert row for every one of the multiple genes
    for r in new_rows_dictionaries:
        df_dictionary = pd.DataFrame([r])
        vus_df = pd.concat([vus_df, df_dictionary], ignore_index=True)

    # sort the dataframe by the 'VUS Id' column in ascending order
    vus_df = vus_df.sort_values(by='VUS Id')

    # reset the indices to maintain a continuous index
    vus_df = vus_df.reset_index(drop=True)

    return vus_df


# create columns for chromosome and chromosome position
def extract_chr_and_position(vus_df: pd.DataFrame):
    # create new columns for chromosome and position and initialize it with empty strings
    vus_df.insert(1, 'Chr', "")
    vus_df.insert(2, 'Position', "")

    # make a copy of the dataframe to be able to iterate through it whilst modifying the original dataframe
    new_vus_df = vus_df.copy()

    # iterate through the dataframe and add entries to the chromosome and position column based on the locus column
    for index, row in new_vus_df.iterrows():
        # get the chromosome and the chromosome position from the locus
        locus = row['Locus']
        locus_arr = locus.split(':')

        vus_df.at[index, 'Chr'] = locus_arr[0].split('chr')[1]
        vus_df.at[index, 'Position'] = locus_arr[1]

    # remove locus column
    vus_df = vus_df.drop(columns=['Locus'])

    return vus_df


def check_for_multiple_genes(vus_df: pd.DataFrame) -> List:
    multiple_genes = []

    for index, row in vus_df.iterrows():
        # get the gene name
        gene_string = row['Gene']
        genes = gene_string.replace(' ', '').split(',')

        # exclude AS1 and LOC from genes
        filtered_genes = [g for g in genes if 'AS1' not in g and 'LOC' not in g]

        # note if variant has multiple genes
        if len(filtered_genes) > 1:
            vus = prep_unprocessed_vus_dict_for_react(row.to_dict())
            multiple_genes.append({'index': index, 'vus': vus, 'genes': filtered_genes})

    return multiple_genes


def extract_sample_ids(sample_ids: str) -> List:
    return re.split(',|;', sample_ids.replace(' ', ''))


def filter_vus(vus_df: pd.DataFrame, one_time_filter_flag: bool) -> InternalResponse:
    if one_time_filter_flag:
        multiple_genes = check_for_multiple_genes(vus_df)
        # TODO: remove other function for multiple genes

        if len(multiple_genes) > 0:
            return InternalResponse({'multiple_genes': multiple_genes}, 200)

    # exclude technical artifacts and CNVs
    vus_df = vus_df[vus_df['Classification'].str.contains('TECHNICAL_ARTIFACT', case=False, regex=True) == False]
    vus_df = vus_df[vus_df['Type'].str.contains('CNV|LONGDEL', case=False, regex=True) == False]

    vus_df = extract_chr_and_position(vus_df)

    return InternalResponse({'multiple_genes': None, 'filtered_df': vus_df}, 200)


def get_rsids(vus_df: pd.DataFrame):
    get_rsids_from_dbsnp_res = get_rsids_from_dbsnp(vus_df)

    if get_rsids_from_dbsnp_res.status != 200:
        current_app.logger.error(
            f'Get RSIDs from dbSNP query failed 500')
        return InternalResponse(None, 500)
    else:
        vus_df = get_rsids_from_dbsnp_res.data

        # vus_df = handle_vus_with_multiple_genes(
        #     vus_df)
        # TODO: is this correct location? should it be included with RSID retrieval process? [needs to be done after rsids]

        # write dataframe to excel file
        # TODO: write to database
        # vus_df.to_excel('rsid_vus.xlsx', index=False)

        return InternalResponse(vus_df, 200)


def get_gene_ids(vus_df: pd.DataFrame) -> pd.DataFrame:
    # create new column for gene id and initialize it with empty strings
    vus_df.insert(3, 'Gene Id', "")

    # make a copy of the dataframe to be able to iterate through it whilst modifying the original dataframe
    new_vus_df = vus_df.copy()

    # iterate through the dataframe
    for index, row in new_vus_df.iterrows():
        # Retrieving gene ids from Gene Annotations table where the VUS is located
        gene_ids: List[GeneAnnotations.gene_id] = db.session.query(GeneAnnotations.id).filter(
            GeneAnnotations.seq_name == row['Chr'],
            row['Position'] >= GeneAnnotations.start_location,
            GeneAnnotations.end_location >= row['Position']
        ).all()

        for gene_id in gene_ids:
            # Retrieving gene name from Gene Attributes table where the VUS is located
            gene_attributes: GeneAttributes = db.session.query(GeneAttributes).filter(
                GeneAttributes.gene_id == gene_id[0],
                GeneAttributes.attribute_name == 'gene_name'
            ).one()

            if gene_attributes.attribute_value not in row['Gene']:
                print(f'The following row:\n {row} \nhas  a gene which does not match the gene found in our database '
                      f'with id:{gene_id} and name:{gene_attributes.attribute_value}')

        if len(gene_ids) > 1:
            print(f'The following row:\n {row} \nhas multiple gene ids: {gene_ids}')
        elif len(gene_ids) == 0:
            print(f'The following row:\n {row} \nhas no gene ids')
        else:
            # TODO: handle multiple genes
            # TODO: handle mismatch genes
            # for now select the first gene
            vus_df.at[index, 'Gene Id'] = gene_ids[0][0]

    return vus_df


def get_external_references(variant_id: int, index: Hashable, vus_df: pd.DataFrame) -> pd.DataFrame:
    # retrieve all external references related to that variant
    external_references: List[ExternalReferences] = db.session.query(ExternalReferences).filter(
        ExternalReferences.variant_id == variant_id
    ).all()

    for ref in external_references:
        if ref.db_type == 'db_snp':
            # retrieve dbsnp entry related to the variant
            dbsnp: DbSnp = db.session.query(DbSnp).filter(
                DbSnp.external_db_snp_id == ref.id
            ).one_or_none()

            if dbsnp.rsid != "NORSID":
                vus_df.at[index, 'RSID'] = dbsnp.rsid
            vus_df.at[index, 'RSID dbSNP verified'] = len(ref.error_msg) == 0
            vus_df.at[index, 'RSID dbSNP errorMsgs'] = ref.error_msg

        elif ref.db_type == 'clinvar':
            # retrieve clinvar entry related to the variant
            clinvar: Clinvar = db.session.query(Clinvar).filter(
                Clinvar.external_clinvar_id == ref.id
            ).one_or_none()

            auto_clinvar_update_id, review_status, classification, last_evaluated = get_last_saved_clinvar_update(clinvar.id)

            # populate the clinvar fields
            vus_df.at[index, 'Clinvar variation id'] = clinvar.variation_id
            vus_df.at[index, 'Clinvar canonical spdi'] = clinvar.canonical_spdi
            vus_df.at[index, 'Clinvar classification'] = classification
            vus_df.at[index, 'Clinvar classification review status'] = review_status
            vus_df.at[index, 'Clinvar classification last eval'] = last_evaluated
            vus_df.at[index, 'Clinvar error msg'] = ref.error_msg

    return vus_df


def check_for_existing_variants(vus_df: pd.DataFrame) -> (pd.DataFrame, pd.DataFrame, List[int]):
    existing_variant_ids = []

    # make a copy of the dataframe to be able to iterate through it whilst modifying the original dataframe
    new_vus_df = vus_df.copy()

    # iterate through the dataframe
    for index, row in new_vus_df.iterrows():
        variant: Variants = get_variant_from_db(row)

        # if variant exists in db
        if variant is not None:
            # populate the remaining fields
            vus_df.at[index, 'Variant Id'] = variant.id
            vus_df.at[index, 'Exists in DB'] = True
            vus_df.at[index, 'Classification'] = variant.classification.value
            # vus_df.at[index, ''] = variant.consequences TODO: add consequence to variant info

            existing_variant_ids.append(variant.id)

            vus_df = get_external_references(variant.id, index, vus_df)

    existing_vus_df = vus_df[vus_df['Exists in DB']]
    vus_df = vus_df[~vus_df['Exists in DB']]

    return existing_vus_df, vus_df, existing_variant_ids


def add_missing_columns(vus_df: pd.DataFrame) -> pd.DataFrame:
    # insert columns for clinvar clinical significance and error messages
    vus_df['Clinvar classification'] = ''
    vus_df['Clinvar classification last eval'] = ""
    vus_df['Clinvar classification review status'] = ""
    vus_df['Clinvar error msg'] = ""
    vus_df['Clinvar canonical spdi'] = ""
    vus_df['Clinvar variation id'] = ""

    # insert columns for dbsnp and error messages
    vus_df['RSID'] = ""
    vus_df['RSID dbSNP verified'] = False
    vus_df['RSID dbSNP errorMsgs'] = ""

    vus_df['Variant Id'] = ""

    # create new column to determine whether the variant has already been stored in the db or not
    vus_df.insert(0, 'Exists in DB', False)

    return vus_df


def get_external_references_for_new_vus(new_vus_df: pd.DataFrame) -> InternalResponse:
    preprocess_and_get_rsids_res = get_rsids(new_vus_df)

    if preprocess_and_get_rsids_res.status != 200:
        current_app.logger.error(
            f'Preprocessing of VUS and retrieval of RSIDs failed 500')
        return InternalResponse({'areRsidsRetrieved:': False, 'isClinvarAccessed': False}, 500)
    else:
        new_vus_df = preprocess_and_get_rsids_res.data

        retrieve_clinvar_variant_classifications_res = retrieve_clinvar_variant_classifications(new_vus_df)

        if retrieve_clinvar_variant_classifications_res.status != 200:
            current_app.logger.error(
                f'Retrieval of ClinVar variant classifications failed 500')
            return InternalResponse({'areRsidsRetrieved:': True, 'isClinvarAccessed': False}, 500)
        else:
            new_vus_df = retrieve_clinvar_variant_classifications_res.data

    return InternalResponse({'areRsidsRetrieved:': True, 'isClinvarAccessed': True,
                             'new_vus_df': new_vus_df}, 200)


def preprocess_vus(vus_df: pd.DataFrame):
    vus_df = add_missing_columns(vus_df)

    if 'Gene Id' not in vus_df.keys():
        vus_df = get_gene_ids(vus_df)

    # check for existing variants
    existing_vus_df, new_vus_df, existing_variant_ids = check_for_existing_variants(vus_df)

    # get external references for new vus
    if len(new_vus_df) > 0:
        get_external_references_for_new_vus_res = get_external_references_for_new_vus(new_vus_df)

        if get_external_references_for_new_vus_res.status != 200:
            current_app.logger.error(
                f'Preprocessing of VUS and retrieval of RSIDs failed 500')
            return InternalResponse(
                {'areRsidsRetrieved:': get_external_references_for_new_vus_res.data['areRsidsRetrieved'],
                 'isClinvarAccessed': get_external_references_for_new_vus_res.data['isClinvarAccessed'],
                 'multiple_genes': None}, 500)
        else:
            new_vus_df = get_external_references_for_new_vus_res.data['new_vus_df']

    # get updated external references for existing vus
    if len(existing_vus_df) > 0:
        get_updated_external_references_for_existing_vus_res = get_updated_external_references_for_existing_vus(existing_vus_df)

        if get_updated_external_references_for_existing_vus_res.status != 200:
            current_app.logger.error(
                f'Preprocessing of VUS and retrieval of RSIDs failed 500')
            return InternalResponse(
                {'areRsidsRetrieved:': True,
                 'isClinvarAccessed': False,
                 'multiple_genes': None}, 500)
        else:
            existing_vus_df = get_updated_external_references_for_existing_vus_res.data['existing_vus_df']

    return InternalResponse({'existing_vus_df': existing_vus_df, 'vus_df': new_vus_df,
                             'existing_variant_ids': existing_variant_ids, 'multiple_genes': None}, 200)


def preprocess_vus_from_file(vus_df: pd.DataFrame, check_for_multiple_genes_flag: bool) -> InternalResponse:
    filter_vus_res = filter_vus(vus_df, check_for_multiple_genes_flag)

    if filter_vus_res.data['multiple_genes'] is not None and len(filter_vus_res.data['multiple_genes']) > 0:
        current_app.logger.info(f'Some VUS contain multiple genes!')
        return InternalResponse({'areRsidsRetrieved:': False, 'isClinvarAccessed': False,
                                 'multiple_genes': filter_vus_res.data['multiple_genes']}, 200)
    else:
        vus_df = filter_vus_res.data['filtered_df']

        return preprocess_vus(vus_df)


def store_new_vus_df_in_db(vus_df: pd.DataFrame) -> List[int]:
    variant_ids = []

    # iterate through the dataframe
    for index, row in vus_df.iterrows():
        # create new variant
        new_variant = Variants(chromosome=row['Chr'], chromosome_position=row['Position'], variant_type=row['Type'],
                               ref=row['Reference'], alt=row['Alt'], classification=Classification.VUS,
                               gene_id=row['Gene Id'], gene_name=row['Gene'])
        # add the new variant to the session
        db.session.add(new_variant)

        db.session.flush()
        variant_ids.append(new_variant.id)

        if len(row['RSID']) > 0 and row['RSID'] != 'NORSID':
            new_dbnsp_external_ref = ExternalReferences(variant_id=new_variant.id,
                                                        db_type='db_snp',
                                                        error_msg=row['RSID dbSNP errorMsgs'])
            db.session.add(new_dbnsp_external_ref)

            db.session.flush()

            new_dbsnp = DbSnp(rsid=row['RSID'],
                              external_db_snp_id=new_dbnsp_external_ref.id)
            db.session.add(new_dbsnp)

        if len(row['Clinvar variation id']) > 0:
            new_clinvar_external_ref = ExternalReferences(variant_id=new_variant.id,
                                                          db_type='clinvar',
                                                          error_msg=row['Clinvar error msg'])
            db.session.add(new_clinvar_external_ref)
            db.session.flush()

            new_clinvar = Clinvar(variation_id=row['Clinvar variation id'],
                                  external_clinvar_id=new_clinvar_external_ref.id,
                                  canonical_spdi=row['Clinvar canonical spdi'])
            db.session.add(new_clinvar)
            db.session.flush()

            store_clinvar_info(new_clinvar.id, row['Clinvar classification'],
                               row['Clinvar classification review status'], row['Clinvar classification last eval'],
                               True)

    return variant_ids


def convert_no_hpo_term_phenotypes_to_array(no_hpo_term_phenotypes_dict: Dict) -> List:
    no_hpo_term_phenotypes_arr = []

    for phenotype in no_hpo_term_phenotypes_dict.keys():
        no_hpo_term_phenotypes_arr.append({'phenotype': phenotype, 'samples': no_hpo_term_phenotypes_dict[phenotype]})

    return no_hpo_term_phenotypes_arr


def update_sample_dict_with_unique_values(sample_id: str, unique_samples_values_dict: Dict, values: List[str]) -> Dict:
    sample_values = unique_samples_values_dict.get(sample_id, [])
    sample_values.extend(values)

    # ensure there are no repeated phenotypes for a given sample
    unique_samples_values_dict[sample_id] = []
    for p in sample_values:
        if p not in unique_samples_values_dict[sample_id]:
            unique_samples_values_dict[sample_id].append(p)

    return unique_samples_values_dict


def convert_dataframe_row_into_array(df_row: pd.Series) -> List:
    array_row = []

    string_row = str(df_row)
    if string_row != 'nan' and len(string_row) > 0:
        array_row = re.split(',', string_row)
        array_row = [p.strip() for p in array_row]

    return array_row


def create_sample_upload_and_sample_entries_in_db(vus_df: pd.DataFrame):
    # retrieve all unique samples and their phenotypes
    unique_samples_phenotypes_dict = {}

    no_hpo_term_phenotypes_dict = {}

    # iterate through the dataframe
    for index, row in vus_df.iterrows():
        # extract sample ids
        sample_ids = extract_sample_ids(str(row['Sample Ids'])) \
 \
        # ontology ids & their names
        phenotype_terms = []

        # extract phenotypes with ids
        if 'Sample Phenotypes With Ids' in vus_df.keys():
            phenotype_terms = row['Sample Phenotypes With Ids']
        # extract phenotypes without ids
        else:
            phenotypes = convert_dataframe_row_into_array(row['Sample Phenotypes'])

            # get phenotype terms (aka ontology ids & their names)
            # retrieve HPO term for each phenotype
            for p in phenotypes:
                hpo_res = get_hpo_term_from_phenotype_name(p)

                if hpo_res.status != 200:
                    samples_no_hpo_terms = no_hpo_term_phenotypes_dict.get(p, [])
                    samples_no_hpo_terms.extend(sample_ids)
                    no_hpo_term_phenotypes_dict[p] = list(set(samples_no_hpo_terms))
                    current_app.logger.error(f'Failed to get HPO term for phenotype {p}')
                else:
                    phenotype_terms.append(hpo_res.data)

        for sample_id in sample_ids:
            # take note of each phenotype for every sample
            unique_samples_phenotypes_dict = update_sample_dict_with_unique_values(sample_id,
                                                                                   unique_samples_phenotypes_dict,
                                                                                   phenotype_terms)

    # store each unique sample in the db & its phenotypes
    for unique_sample_id in unique_samples_phenotypes_dict.keys():
        # check if sample already exists
        sample: Samples = db.session.query(Samples).filter(Samples.id == unique_sample_id).one_or_none()

        # if the sample is new, add it to database
        if sample is None:
            sample = Samples(id=unique_sample_id, genome_version='GRCh37')
            db.session.add(sample)

            db.session.flush()

        sample_ontology_term_ids = [o.ontology_term_id for o in sample.ontology_term]

        # append phenotypes to the respective sample, if the sample does not already have that phenotype
        for phenotype_term in unique_samples_phenotypes_dict[sample.id]:
            if phenotype_term['ontologyId'] not in sample_ontology_term_ids:
                append_phenotype_to_sample(sample, phenotype_term)

    db.session.flush()

    no_hpo_term_phenotypes = convert_no_hpo_term_phenotypes_to_array(no_hpo_term_phenotypes_dict)

    return InternalResponse({'no_hpo_term_phenotypes': no_hpo_term_phenotypes}, 200)


def store_acmg_rules_for_variant(are_rules_with_ids: bool, vus_df: pd.DataFrame, variant_ids: List[int]):
    # iterate through the dataframe
    for index, row in vus_df.iterrows():
        variant_id = variant_ids[int(index)]

        # extract acmg rules with ids
        if are_rules_with_ids:
            acmg_rules = row['ACMG Rules With Ids']
        # extract acmg rules without ids
        else:
            acmg_rules = []
            acmg_rules_input = list(set(convert_dataframe_row_into_array(row['ACMG Rules'])))
            for rule in acmg_rules_input:
                acmg_rule: AcmgRules = db.session.query(AcmgRules).filter(AcmgRules.rule_name == rule).first()
                acmg_rules.append({'id': acmg_rule.id, 'name': acmg_rule.rule_name})
            # TODO include mapping of acmg rules

        # store acmg rules to the respective variant, if the variant does not already have that acmg rule
        new_added_acmg_rule_ids = []
        for rule in acmg_rules:
            # check if variant already has this acmg rule
            variants_acmg_rule: VariantsAcmgRules | None = db.session.query(
                VariantsAcmgRules).filter(VariantsAcmgRules.acmg_rule_id == rule['id'],
                                          VariantsAcmgRules.variant_id == variant_id).one_or_none()

            if variants_acmg_rule is None:
                # store acmg rules for the given sample
                new_variants_acmg_rule = VariantsAcmgRules(variant_id=variant_id,
                                                           acmg_rule_id=rule['id'],
                                                           rule_name=rule['name'])
                db.session.add(new_variants_acmg_rule)
                new_added_acmg_rule_ids.append(new_variants_acmg_rule.acmg_rule_id)

        # ensure classification exists
        if row['Classification'] in Classification.__members__:
            classification = row['Classification']
        else:
            classification = 'VUS'

        # create review to record variant creation
        new_review: Reviews = Reviews(variant_id=variant_id, scientific_member_id=current_user.id,
                                      date_added=datetime.now(),
                                      classification=classification, classification_reason=None)
        new_review.acmg_rules = db.session.query(AcmgRules).filter(AcmgRules.id.in_(new_added_acmg_rule_ids)).all()
        db.session.add(new_review)


def store_upload_details_for_variant_sample(file_upload: FileUploads | None, is_file_upload: bool, sample_id: str,
                                            variant_id: int):
    # create sample upload entry
    if is_file_upload:
        upload_type = 'file'
    else:
        upload_type = 'manual'

    new_variants_sample_upload = VariantsSamplesUploads(date_uploaded=datetime.now(),
                                                        scientific_member_id=current_user.id,
                                                        upload_type=upload_type, sample_id=sample_id,
                                                        variant_id=variant_id)

    if is_file_upload and file_upload is not None:
        new_variants_sample_upload.file_upload = file_upload

    db.session.add(new_variants_sample_upload)
    db.session.flush()

    if not is_file_upload:
        # create manual uploads entry
        new_manual_upload = ManualUploads(variants_samples_uploads_manual_id=new_variants_sample_upload.id)
        db.session.add(new_manual_upload)


def store_variant_sample_relations_in_db(vus_df: pd.DataFrame, variant_ids: List[int], file: FileStorage | None,
                                         is_file_upload: bool):
    if is_file_upload:
        # create entry for file in the db
        new_file_upload = FileUploads(filename=file.filename)
        db.session.add(new_file_upload)

        db.session.flush()
    else:
        new_file_upload = None

    # iterate through the dataframe
    for index, row in vus_df.iterrows():
        variant_id = variant_ids[int(index)]

        alt = row['Alt']
        vus_genotype = row['Genotype']
        hgvs = row['HGVS']

        genotype_split = vus_genotype.split('/')

        genotype = Genotype.HETEROZYGOUS

        if genotype_split[0] == alt and genotype_split[1] == alt:
            genotype = Genotype.HOMOZYGOUS

        samples = extract_sample_ids(str(row['Sample Ids']))

        for sample_id in samples:
            # check if variants_samples entry already exists
            existing_variants_samples: VariantsSamples = db.session.query(VariantsSamples).filter(
                VariantsSamples.variant_id == variant_id, VariantsSamples.sample_id == sample_id).one_or_none()

            # if variants_samples entry does not exist, add new entry
            if existing_variants_samples is None:
                # check if HGVS exists
                variant_hgvs: VariantHgvs = db.session.query(VariantHgvs).filter(VariantHgvs.variant_id == variant_id, VariantHgvs.hgvs == hgvs).one_or_none()

                if variant_hgvs is None:
                    # create new HGVS entry
                    variant_hgvs = VariantHgvs(variant_id=variant_id, hgvs=hgvs, is_updated=False)
                    db.session.add(variant_hgvs)
                    db.session.flush()

                new_variants_samples = VariantsSamples(variant_id=variant_id, sample_id=sample_id, variant_hgvs_id=variant_hgvs.id ,genotype=genotype)
                db.session.add(new_variants_samples)

            # store the upload details related to this variant & sample
            store_upload_details_for_variant_sample(new_file_upload, is_file_upload, sample_id, variant_id)


def store_vus_info_in_db(existing_vus_df: pd.DataFrame, existing_variant_ids: List[str], new_vus_df: pd.DataFrame,
                         file: FileStorage | None, is_file_upload: bool) -> Response:
    # join existing vus_df with the new vus
    all_vus_df = pd.concat([existing_vus_df, new_vus_df], axis=0)

    # store file and samples entries
    create_sample_upload_and_sample_entries_in_db_res = create_sample_upload_and_sample_entries_in_db(all_vus_df)

    # store phenotypes (and respective samples) which do not have an exact match to an HPO term
    no_hpo_term_phenotypes = create_sample_upload_and_sample_entries_in_db_res.data['no_hpo_term_phenotypes']

    # store those vus that do not already exist in the db
    variant_ids = store_new_vus_df_in_db(new_vus_df)
    new_vus_df['Variant Id'] = variant_ids

    # override to include newly added variants' ids
    all_vus_df = pd.concat([existing_vus_df, new_vus_df], axis=0)

    # retrieve and store the publications (user links & LitVar) of new variants
    retrieve_and_store_variant_pub_res = retrieve_and_store_variant_publications(new_vus_df, False)

    if retrieve_and_store_variant_pub_res.status != 200:
        current_app.logger.error(
            f'Get publications for new variants failed 500')

    # retrieve and store the publications (user links & LitVar) of existing variants
    retrieve_and_store_existing_variant_pub_res = retrieve_and_store_variant_publications(existing_vus_df, True)

    if retrieve_and_store_existing_variant_pub_res.status != 200:
        current_app.logger.error(
            f'Get publications for existing variants failed 500')

    # join existing vus_df's variant ids with the new vus variant ids
    all_variant_ids = existing_variant_ids + variant_ids

    # store the acmg rules related to this variant
    store_acmg_rules_for_variant('ACMG Rules With Ids' in all_vus_df.keys(), all_vus_df, all_variant_ids)

    # store variants-samples entries in db
    store_variant_sample_relations_in_db(all_vus_df, all_variant_ids, file, is_file_upload)

    try:
        # Commit the session to persist changes to the database
        db.session.commit()
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since insertion of entries in DB failed due to error: {e}')

        response = {'isSuccess': False, 'areRsidsRetrieved:': True, 'isClinvarAccessed': True}
        return Response(json.dumps(response), 200)

    # update column names to camelCase format
    new_vus_df = prep_vus_df_for_react(all_vus_df)

    # convert df to list
    vus_list = convert_df_to_list(new_vus_df)

    return Response(
        json.dumps({'isSuccess': True, 'areRsidsRetrieved:': True, 'isClinvarAccessed': True,
                    'vusList': vus_list, 'noHpoTermPhenotypes': no_hpo_term_phenotypes}), 200)


def handle_vus_file(file: FileStorage, multiple_genes_selection: List) -> Response:
    vus_df = pd.read_excel(file, header=0)  # TODO: check re header
    current_app.logger.info(f'Number of VUS found in file: {len(vus_df)}')

    # handle multiple genes selection
    if len(multiple_genes_selection) > 0:
        for selection in multiple_genes_selection:
            vus_df.at[int(selection['index']), 'Gene'] = selection['gene']

    # flag indicating that check is required for multiple genes for a single variant
    one_time_filter_flag = len(multiple_genes_selection) == 0

    preprocess_vus_res = preprocess_vus_from_file(vus_df, one_time_filter_flag)

    if preprocess_vus_res.status != 200:
        response = preprocess_vus_res.data
        response['isSuccess'] = False
        return Response(json.dumps(response), 200)
    elif (one_time_filter_flag and preprocess_vus_res.data['multiple_genes'] is not None
          and len(preprocess_vus_res.data['multiple_genes']) > 0):
        return Response(json.dumps({'isSuccess': True, 'multipleGenes': preprocess_vus_res.data['multiple_genes']}),
                        200, mimetype='application/json')
    else:
        existing_vus_df = preprocess_vus_res.data['existing_vus_df']
        existing_variant_ids = preprocess_vus_res.data['existing_variant_ids']

        vus_df = preprocess_vus_res.data['vus_df']

        return store_vus_info_in_db(existing_vus_df, existing_variant_ids, vus_df, file, True)


def handle_vus_from_form(vus_df: pd.DataFrame) -> Response:
    # adjust existing column names - rename the existing DataFrame
    vus_df.rename(columns={'chromosome': 'Chr', 'chromosomePosition': 'Position', 'type': 'Type',
                           'refAllele': 'Reference', 'altAllele': 'Alt', 'classification': Classification.VUS,
                           'gene': 'Gene', 'geneId': 'Gene Id', 'genotype': 'Genotype', 'samples': 'Sample Ids',
                           'phenotypes': 'Sample Phenotypes With Ids', 'acmgRules': 'ACMG Rules With Ids',
                           'hgvs': 'HGVS'},
                  inplace=True)

    preprocess_vus_res = preprocess_vus(vus_df)

    if preprocess_vus_res.status != 200:
        response = preprocess_vus_res.data
        response['isSuccess'] = False
        return Response(json.dumps(response), 200)
    else:
        existing_vus_df = preprocess_vus_res.data['existing_vus_df']
        existing_variant_ids = preprocess_vus_res.data['existing_variant_ids']

        vus_df = preprocess_vus_res.data['vus_df']
        vus_df['Classification'] = 'VUS'

        return store_vus_info_in_db(existing_vus_df, existing_variant_ids, vus_df, None, False)
