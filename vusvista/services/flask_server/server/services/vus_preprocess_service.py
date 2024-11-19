import html
import math
import traceback
from datetime import datetime
from io import BytesIO
from typing import List, Hashable, Dict

from flask import current_app, Response
import pandas as pd
from Bio import Entrez
from flask_login import current_user
from sqlalchemy import and_
from sqlalchemy.exc import SQLAlchemyError
import json
import re

from server import db, file_upload_tasks
from server.helpers.data_helper import prep_vus_df_for_react, convert_df_to_list, prep_unprocessed_vus_dict_for_react
from server.helpers.db_access_helper import get_variant_from_db
from server.models import Variants, GeneAnnotations, GeneAttributes, DbSnp, \
    Clinvar, ExternalReferences, FileUploads, Genotype, \
    AcmgRules, VariantsAcmgRules, Classification, Reviews, FileUploadEvents
from server.responses.internal_response import InternalResponse
from server.services.consequence_service import get_consequences_for_new_vus
from server.services.dbsnp_service import get_rsids_from_dbsnp

from server.services.clinvar_service import retrieve_clinvar_variant_classifications, get_updated_external_references_for_existing_vus, store_clinvar_info
from server.services.phenotype_service import get_hpo_term_from_phenotype_name
from server.services.samples_service import add_new_sample_to_db
from server.services.variants_samples_service import store_upload_details_for_variant_sample, add_variant_sample_to_db
from server.services.view_vus_service import get_last_saved_clinvar_update, retrieve_vus_summaries_from_db
from server.services.publications_service import retrieve_and_store_variant_publications

Entrez.email = "esther.spiteri.18@um.edu.mt"

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

        # to handle e.g. chr10:88439188 A⇒G
        vus_df.at[index, 'Position'] = locus_arr[1].split(' ')[0].split('_')[0]

    # remove locus column
    vus_df = vus_df.drop(columns=['Locus'])

    return vus_df


def get_filtered_genes(gene_row: str):
    # get the gene name
    gene_string = gene_row
    genes = gene_string.replace(' ', '').split(',')

    # exclude AS1 and LOC from genes
    filtered_genes = [g for g in genes if 'AS1' not in g and 'LOC' not in g]

    return filtered_genes


def check_for_multiple_genes(vus_df: pd.DataFrame) -> List:
    multiple_genes = []

    # make a copy of the dataframe to be able to iterate through it whilst modifying the original dataframe
    vus_df_copy = vus_df.copy()

    try:
        for index, row in vus_df_copy.iterrows():
            # exclude AS1 and LOC from genes
            filtered_genes = get_filtered_genes(row['Gene'])

            # note if variant has multiple genes
            if len(filtered_genes) > 1:
                vus = prep_unprocessed_vus_dict_for_react(row.to_dict())
                multiple_genes.append({'index': index, 'vus': vus, 'genes': filtered_genes})
    except Exception as e:
        current_app.logger.error(current_app.logger.error(traceback.format_exc()))

    return multiple_genes


def check_for_existing_genes(vus_df: pd.DataFrame) -> List:
    genes_not_found_in_db = []

    for index, row in vus_df.iterrows():
        # exclude AS1 and LOC from genes
        filtered_genes = get_filtered_genes(row['Gene'])

        # assume that only one gene is left after filtering (& given that multiple gene check has been executed)
        gene_db: GeneAttributes = db.session.query(GeneAttributes).filter(and_(GeneAttributes.attribute_name.in_(['gene_name', 'gene_id']), GeneAttributes.attribute_value == filtered_genes[0])).one_or_none()

        if gene_db is None:
            vus = prep_unprocessed_vus_dict_for_react(row.to_dict())
            genes_not_found_in_db.append({"index": index, "vus": vus, "gene": row['Gene']})

    return genes_not_found_in_db


def extract_sample_ids(sample_ids: str) -> List:
    return re.split(',|;', sample_ids.replace(' ', ''))


def filter_vus(vus_df: pd.DataFrame) -> InternalResponse:
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

        return InternalResponse(vus_df, 200)


def get_gene_ids(vus_df: pd.DataFrame) -> pd.DataFrame:
    # create new column for gene id and initialize it with empty strings
    vus_df.insert(3, 'Gene Id', "")

    # make a copy of the dataframe to be able to iterate through it whilst modifying the original dataframe
    new_vus_df = vus_df.copy()

    # iterate through the dataframe
    for index, row in new_vus_df.iterrows():
        filtered_genes = [g.upper() for g in get_filtered_genes(row['Gene'])]
        is_matching_gene_found = False

        # Retrieving gene ids from Gene Annotations table where the VUS is located
        gene_ids: List[GeneAnnotations.gene_id] = db.session.query(GeneAnnotations.id).filter(
            GeneAnnotations.seq_name == row['Chr'],
            row['Position'] >= GeneAnnotations.start_location,
            GeneAnnotations.end_location >= row['Position']
        ).all()

        for gene_id in gene_ids:
            # Retrieving gene name from Gene Attributes table where the VUS is located
            gene_attributes: List[GeneAttributes] = db.session.query(GeneAttributes).filter(
                GeneAttributes.gene_id == gene_id[0],
                (GeneAttributes.attribute_name == 'gene_name') |
                (GeneAttributes.attribute_name == 'gene_id')
            ).all()

            for a in gene_attributes:
                if a.attribute_value in filtered_genes:
                    is_matching_gene_found = True
                    vus_df.at[index, 'Gene Id'] = gene_ids[0][0]
                    vus_df.at[index, 'Gene'] = a.attribute_value
                    break

        # ignore variant location and get gene id of inputted gene
        if not is_matching_gene_found:
            current_app.logger.error(f'The following row:\n {row} \nhas a gene which does not match the genes found in our database '
                  f'with ids:{gene_ids} that match to the variant\'s coordinates. Checking for gene\'s with an exact gene name match.')
            gene_attributes: List[GeneAttributes] = db.session.query(GeneAttributes).filter(
                (and_(GeneAttributes.attribute_name == 'gene_name', GeneAttributes.attribute_value.in_(filtered_genes))) |
                (and_(GeneAttributes.attribute_name == 'gene_id', GeneAttributes.attribute_value.in_(filtered_genes)))
            ).all()

            for a in gene_attributes:
                if a.attribute_value in filtered_genes:
                    vus_df.at[index, 'Gene Id'] = a.gene_id
                    vus_df.at[index, 'Gene'] = a.attribute_value
                    break

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

    vus_df['Consequence'] = ""

    # create new column to determine whether the variant has already been stored in the db or not
    vus_df.insert(0, 'Exists in DB', False)

    return vus_df


def get_external_references_for_new_vus(new_vus_df: pd.DataFrame) -> InternalResponse:
    try:
        preprocess_and_get_rsids_res = get_rsids(new_vus_df)
    except Exception as e:
        current_app.logger.error(f'RSID retrieval failed because of the error: {e}')
        return InternalResponse({'areRsidsRetrieved': False, 'isClinvarAccessed': False}, 500)

    if preprocess_and_get_rsids_res.status != 200:
        current_app.logger.error(
            f'Preprocessing of VUS and retrieval of RSIDs failed when getting external ref for new vus 500')
        return InternalResponse({'areRsidsRetrieved': False, 'isClinvarAccessed': False}, 500)
    else:
        new_vus_df = preprocess_and_get_rsids_res.data

        retrieve_clinvar_variant_classifications_res = retrieve_clinvar_variant_classifications(new_vus_df)

        if retrieve_clinvar_variant_classifications_res.status != 200:
            current_app.logger.error(
                f'Retrieval of ClinVar variant classifications failed 500')
            return InternalResponse({'areRsidsRetrieved': True, 'isClinvarAccessed': False}, 500)
        else:
            new_vus_df = retrieve_clinvar_variant_classifications_res.data

    return InternalResponse({'areRsidsRetrieved': True, 'isClinvarAccessed': True,
                             'new_vus_df': new_vus_df}, 200)


def get_consequences_for_new_uploads(vus_df: pd.DataFrame) -> pd.DataFrame:
    hgvs = []
    for var in vus_df.iterrows():
        if not (isinstance(var[1]["HGVS"], float) and math.isnan(var[1]["HGVS"])) and var[1]['HGVS'] is not None:
            hgvs.append(var[1]["HGVS"])

    if len(hgvs) > 0:
        # get variant consequences through HGVS
        get_consequences_for_new_vus_res = get_consequences_for_new_vus(hgvs)

        if get_consequences_for_new_vus_res.status != 200:
            current_app.logger.error(
                f'Preprocessing of VUS and retrieval of variant consequences failed 500')
        else:
            hgvs_dict = get_consequences_for_new_vus_res.data['consequences_dict']

            # make a copy of the dataframe to be able to iterate through it whilst modifying the original dataframe
            vus_df_copy = vus_df.copy()
            for index, row in vus_df_copy.iterrows():
                vus_df.at[index, 'Consequence'] = hgvs_dict.get(row['HGVS'], "")

    return vus_df


def preprocess_vus(vus_df: pd.DataFrame):
    vus_df = add_missing_columns(vus_df)

    if 'Gene Id' not in vus_df.keys():
        vus_df = get_gene_ids(vus_df)

    # check for existing variants
    existing_vus_df, new_vus_df, existing_variant_ids = check_for_existing_variants(vus_df)

    # handle new vus
    if len(new_vus_df) > 0:
        # get external references
        get_external_references_for_new_vus_res = get_external_references_for_new_vus(new_vus_df)

        if get_external_references_for_new_vus_res.status != 200:
            current_app.logger.error(
                f'Preprocessing of VUS and retrieval of RSIDs failed when preprocessing vus 500')
            return InternalResponse(
                {'areRsidsRetrieved:': get_external_references_for_new_vus_res.data['areRsidsRetrieved'],
                 'isClinvarAccessed': get_external_references_for_new_vus_res.data['isClinvarAccessed'],
                 'multiple_genes': None}, 500)
        else:
            new_vus_df = get_external_references_for_new_vus_res.data['new_vus_df']

            new_vus_df = get_consequences_for_new_uploads(new_vus_df)

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

            existing_vus_df = get_consequences_for_new_uploads(existing_vus_df)

    return InternalResponse({'existing_vus_df': existing_vus_df, 'vus_df': new_vus_df,
                             'existing_variant_ids': existing_variant_ids, 'multiple_genes': None}, 200)


def preprocess_vus_from_file(vus_df: pd.DataFrame) -> InternalResponse:
    filter_vus_res = filter_vus(vus_df)

    vus_df = filter_vus_res.data['filtered_df']

    return preprocess_vus(vus_df)


def store_new_vus_df_in_db(vus_df: pd.DataFrame) -> List[int]:
    variant_ids = []

    # iterate through the dataframe
    for index, row in vus_df.iterrows():
        variant = get_variant_from_db(row)

        if variant is None:
            # create new variant
            new_variant = Variants(chromosome=row['Chr'], chromosome_position=row['Position'], variant_type=html.unescape(row['Type']),
                                   ref=row['Reference'], alt=row['Alt'], classification=Classification.VUS,
                                   gene_id=row['Gene Id'], gene_name=row['Gene'])
            # add the new variant to the session
            db.session.add(new_variant)

            db.session.flush()
            variant_ids.append(new_variant.id)

            if str(row['RSID']) != 'nan' and len(row['RSID']) > 0 and row['RSID'] != 'NORSID':
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
        else:
            variant_ids.append(variant.id)
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


def convert_phenotypes_dataframe_row_into_array(phenotypes: pd.Series) -> List:
    phenotypes_arr = []

    phenotypes_str = str(phenotypes)
    if phenotypes_str != 'nan' and len(phenotypes_str) > 0:
        phenotypes_arr = re.split(',', phenotypes_str)
        phenotypes_arr = [p.strip() for p in phenotypes_arr]

    return phenotypes_arr


def convert_acmg_rules_dataframe_row_into_array(acmg_rules: pd.Series) -> List:
    acmg_rules_arr = []

    acmg_rules_str = str(acmg_rules).replace(' ', '')
    if acmg_rules_str != 'nan' and len(acmg_rules_str) > 0:
        acmg_rules_arr = re.split(',', acmg_rules_str)

        acmg_rules_arr = [r.split('_')[0] for r in acmg_rules_arr]

    return acmg_rules_arr


def create_sample_upload_and_sample_entries_in_db(vus_df: pd.DataFrame):
    # retrieve all unique samples and their phenotypes
    unique_samples_phenotypes_dict = {}

    no_hpo_term_phenotypes_dict = {}

    # iterate through the dataframe
    for index, row in vus_df.iterrows():
        # extract unique sample ids
        sample_ids = list(set(extract_sample_ids(str(row['Sample Ids']))))

        # ontology ids & their names
        phenotype_terms = []

        # extract phenotypes with ids
        if 'Sample Phenotypes With Ids' in vus_df.keys():
            phenotype_terms = row['Sample Phenotypes With Ids']
        # extract phenotypes without ids
        else:
            phenotypes = convert_phenotypes_dataframe_row_into_array(row['Sample Phenotypes'])

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
        add_new_sample_to_db(unique_sample_id, unique_samples_phenotypes_dict[unique_sample_id])

    db.session.flush()

    no_hpo_term_phenotypes = convert_no_hpo_term_phenotypes_to_array(no_hpo_term_phenotypes_dict)

    return InternalResponse({'no_hpo_term_phenotypes': no_hpo_term_phenotypes}, 200)


def store_acmg_rules_for_variant(are_rules_with_ids: bool, vus_df: pd.DataFrame, variant_ids: List[int], user_id: int):
    # iterate through the dataframe
    for index, row in vus_df.iterrows():
        variant_id = variant_ids[int(index)]

        # extract acmg rules with ids
        if are_rules_with_ids:
            acmg_rules = row['ACMG Rules With Ids']
        # extract acmg rules without ids
        else:
            acmg_rules = []
            acmg_rules_input = list(set(convert_acmg_rules_dataframe_row_into_array(row['ACMG Rules'])))
            for rule in acmg_rules_input:
                acmg_rule: AcmgRules = db.session.query(AcmgRules).filter(AcmgRules.rule_name == rule).first()
                acmg_rules.append({'id': acmg_rule.id, 'name': acmg_rule.rule_name})

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
        new_review: Reviews = Reviews(variant_id=variant_id, scientific_member_id=user_id,
                                      date_added=datetime.now(),
                                      classification=classification, classification_reason=None)
        new_review.acmg_rules = db.session.query(AcmgRules).filter(AcmgRules.id.in_(new_added_acmg_rule_ids)).all()
        db.session.add(new_review)


def store_variant_sample_relations_in_db(vus_df: pd.DataFrame, variant_ids: List[int], filename: str | None,
                                         is_file_upload: bool, user_id: int):
    if is_file_upload:
        # create entry for file in the db
        new_file_upload = FileUploads(filename=filename)
        db.session.add(new_file_upload)

        db.session.flush()
    else:
        new_file_upload = None

    # iterate through the dataframe
    for index, row in vus_df.iterrows():
        variant_id = variant_ids[int(index)]

        vus_genotype = row['Genotype']
        if '/' in vus_genotype:
            genotype_split = vus_genotype.split('/')
            ref = genotype_split[0]
            alt = genotype_split[1]

            if ref == alt:
                genotype = Genotype.HOMOZYGOUS
            else:
                genotype = Genotype.HETEROZYGOUS
        else:
            genotype = Genotype.HETEROZYGOUS

            if 'homozygous' in vus_genotype.lower():
                genotype = Genotype.HOMOZYGOUS

        samples = list(set(extract_sample_ids(str(row['Sample Ids']))))

        hgvs = row['HGVS']

        for sample_id in samples:
            # store the variants sample
            add_variant_sample_to_db(variant_id, sample_id, hgvs, genotype.value, row['Consequence'])

            # store the upload details related to this variant & sample
            store_upload_details_for_variant_sample(new_file_upload, is_file_upload, sample_id, variant_id, user_id)


def store_vus_info_in_db(existing_vus_df: pd.DataFrame, existing_variant_ids: List[str], new_vus_df: pd.DataFrame,
                         filename: str | None, is_file_upload: bool, user_id: int) -> InternalResponse:
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

    # retrieve and store the publications (user links & LitVar) of all variants
    retrieve_and_store_variant_pub_res = retrieve_and_store_variant_publications(all_vus_df)

    if retrieve_and_store_variant_pub_res.status != 200:
        current_app.logger.error(
            f'Get publications for variants failed 500')

    # join existing vus_df's variant ids with the new vus variant ids
    all_variant_ids = existing_variant_ids + variant_ids

    # store the acmg rules related to this variant
    store_acmg_rules_for_variant('ACMG Rules With Ids' in all_vus_df.keys(), all_vus_df, all_variant_ids, user_id)

    # store variants-samples entries in db
    store_variant_sample_relations_in_db(all_vus_df, all_variant_ids, filename, is_file_upload, user_id)

    try:
        # Commit the session to persist changes to the database
        db.session.commit()
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since insertion of entries in DB failed due to error: {e}')
        return InternalResponse({'isSuccess': False, 'areRsidsRetrieved:': True, 'isClinvarAccessed': True}, 500)

    # update column names to camelCase format
    new_vus_df = prep_vus_df_for_react(all_vus_df)

    # convert df to list
    vus_list = convert_df_to_list(new_vus_df)
    return InternalResponse({'isSuccess': True, 'areRsidsRetrieved:': True, 'isClinvarAccessed': True,
                    'vusList': vus_list, 'noHpoTermPhenotypes': no_hpo_term_phenotypes}, 200)


def handle_vus_file(task_id: int, vus_df: pd.DataFrame, filename: str, user_id: int):
    try:
        vus_df_copy = vus_df.copy()

        # handle empty ref & alt alleles
        for index, row in vus_df_copy.iterrows():
            ref = row['Reference']
            if isinstance(ref, float) and math.isnan(ref):
                vus_df.at[index, 'Reference'] = None

            alt = row['Alt']
            if isinstance(alt, float) and math.isnan(alt):
                vus_df.at[index, 'Alt'] = None

        preprocess_vus_res = preprocess_vus_from_file(vus_df)

        if preprocess_vus_res.status != 200:
            file_upload_tasks[task_id] = {'taskId': task_id, 'isSuccess': False, 'filename': filename}
        else:
            existing_vus_df = preprocess_vus_res.data['existing_vus_df']
            existing_variant_ids = preprocess_vus_res.data['existing_variant_ids']

            vus_df = preprocess_vus_res.data['vus_df']

            store_vus_info_in_db_res = store_vus_info_in_db(existing_vus_df, existing_variant_ids, vus_df, filename, True, user_id)

            # get variant summaries
            variant_ids = [v['id'] for v in store_vus_info_in_db_res.data['vusList']]

            variants: List[Variants] = db.session.query(Variants).filter(Variants.id.in_(variant_ids)).all()
            vus_list = retrieve_vus_summaries_from_db(variants)

            task_status = store_vus_info_in_db_res.data
            task_status['vusList'] = vus_list
            task_status['existingVariantIds'] = existing_variant_ids

            file_upload_tasks[task_id] = task_status
            file_upload_tasks[task_id]['filename'] = filename
            file_upload_tasks[task_id]['taskId'] = task_id
    except Exception as error:
        current_app.logger.error(traceback.format_exc())
        current_app.logger.error(error)
        file_upload_tasks[task_id] = {'taskId': task_id, 'isSuccess': False, 'filename': filename}


def handle_vus_from_form(vus_df: pd.DataFrame) -> Response:
    # adjust existing column names - rename the existing DataFrame
    vus_df.rename(columns={'chromosome': 'Chr', 'chromosomePosition': 'Position', 'type': 'Type',
                           'refAllele': 'Reference', 'altAllele': 'Alt', 'classification': Classification.VUS,
                           'gene': 'Gene', 'geneId': 'Gene Id', 'genotype': 'Genotype', 'samples': 'Sample Ids',
                           'phenotypes': 'Sample Phenotypes With Ids', 'acmgRules': 'ACMG Rules With Ids',
                           'hgvs': 'HGVS', 'rsid': 'RSID_', 'literatureLinks': 'Literature Links'},
                  inplace=True)

    vus_df_copy = vus_df.copy()

    # handle empty ref & alt alleles
    for index, row in vus_df_copy.iterrows():
        ref = row['Reference']
        if isinstance(ref, float) and math.isnan(ref):
            vus_df.at[index, 'Reference'] = None

        alt = row['Alt']
        if isinstance(alt, float) and math.isnan(alt):
            vus_df.at[index, 'Alt'] = None

        if row['HGVS'] == "":
            vus_df.at[index, 'HGVS'] = None

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

        store_vus_info_in_db_res = store_vus_info_in_db(existing_vus_df, existing_variant_ids, vus_df, None, False, current_user.id)

        return Response(json.dumps(store_vus_info_in_db_res.data), 200)


def scheduled_file_upload_events():
    file_upload_events: List[FileUploadEvents] = db.session.query(FileUploadEvents).filter(FileUploadEvents.date_processed.is_(None)).all()

    if file_upload_events:
        for e in file_upload_events:
            current_app.logger.info(f'Handling file Upload with Id {e.id}')

            # Use BytesIO to read the binary content
            byte_stream = BytesIO(e.file_data)

            # Read the Excel content into a DataFrame
            vus_df = pd.read_excel(byte_stream, engine='openpyxl')

            handle_vus_file(e.id, vus_df, e.file_name, e.scientific_member_id)

            e.date_processed = datetime.now()
        try:
            # Commit the session to persist changes to the database
            db.session.commit()
        except SQLAlchemyError as e:
            # Changes were rolled back due to an error
            db.session.rollback()

            current_app.logger.error(
                f'Rollback carried out since storing updated file upload events for ids {[e.id for e in file_upload_events]} in DB failed due to error: {e}')
