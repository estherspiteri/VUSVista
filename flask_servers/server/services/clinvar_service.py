import os
from datetime import datetime

import requests
import xmltodict
from flask import current_app, Response, render_template
import pandas as pd
import time
import json
from Bio import Entrez
from flask_mail import Message
from requests import RequestException

from typing import Dict, List, Tuple

from sqlalchemy import desc
from sqlalchemy.exc import SQLAlchemyError

from server import db
from server.config import mail
from server.helpers.data_helper import get_variant_summary
from server.models import Variants, ExternalReferences, Clinvar, AutoClinvarEvalDates, AutoClinvarUpdates, \
    ScientificMembers
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
def retrieve_clinvar_dict(clinvar_variation_id: str):
    # retrieve clinvar info for a single variant
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id={clinvar_variation_id}&from_esearch=true"

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


def retrieve_multiple_clinvar_dict(clinvar_ids: List[str]):  # TODO: use this to replace retrieve_clinvar_dict
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
def compare_clinvar_variant_with_expected_variant(genome_version: str, classified_record_dict, gene: str,
                                                  chr: str, chr_pos: str) -> tuple[bool, str]:
    # only compare if Clinvar entry has classification (variationId: 478905)
    if classified_record_dict is None:
        return False, (f"Unable to verify this Clinvar entry!")
    else:
        clinvar_allele = classified_record_dict.get('SimpleAllele')

        if clinvar_allele:
            # compare expected and retrieved variant's gene
            clinvar_genes = clinvar_allele.get('GeneList')

            # if none of the genes match the expected gene
            if isinstance(clinvar_genes['Gene'], list):
                if gene not in [g['@Symbol'] for g in clinvar_genes['Gene']]:
                    return False, (f"None of the gene names {[g['@Symbol'] for g in clinvar_genes['Gene']]} "
                                   f"match the expected gene name {gene}!")
            else:
                if gene != clinvar_genes['Gene']['@Symbol']:
                    return False, (f"None of the gene names {clinvar_genes['Gene']['@Symbol']} "
                                   f"match the expected gene name {gene}!")

            clinvar_locations = clinvar_allele.get('Location').get('SequenceLocation')
            for loc in clinvar_locations:
                if loc.get('@Assembly') == genome_version:
                    # compare expected and retrieved variant's chromosome
                    if loc.get('@Chr') != chr:
                        return False, f"Chromosome {loc.get('@Chr')} does not match the expected chromosome {chr}!"
                    # compare expected and retrieved variant's chromosome position
                    elif loc.get('@start') != str(chr_pos):
                        return False, (
                            f"Chromosome start position {loc.get('@start')} does not match the expected chromosome "
                            f"position {chr_pos}!")
                    break
                # TODO compare ref & alt alleles
        else:
            return False, (f"Unable to verify this Clinvar entry!")

    return True, ''


# Retrieve the ClinVar variant's clinical significance.
def extract_clinvar_germline_classification(classified_record_dict: Dict):
    # only compare if Clinvar entry has classification (variationId: 478905)
    if classified_record_dict is None:
        return {'description': None, 'last_evaluated': None, 'review_status': None}
    else:
        germline_classification = classified_record_dict.get('Classifications').get('GermlineClassification')

        last_eval = ""
        if germline_classification.get('@DateLastEvaluated') is not None and germline_classification.get(
                '@DateLastEvaluated') != "1/01/01 00:00":
            last_eval = germline_classification.get('@DateLastEvaluated').replace('-', '/') + ' 00:00'

        return {'description': germline_classification.get('Description'), 'last_evaluated': last_eval,
                'review_status': germline_classification.get('ReviewStatus')}


# Retrieve the ClinVar variant's canonical SPDI.
def extract_clinvar_canonical_spdi(classified_record_dict: Dict):
    if classified_record_dict is None or classified_record_dict.get('SimpleAllele') is None:
        return None

    canonical_spdi = (classified_record_dict
                      .get('SimpleAllele').get('CanonicalSPDI'))

    return canonical_spdi


# Retrieve the ClinVar variation's id.
def extract_clinvar_variation_id(clinvar_dict: Dict):
    return clinvar_dict.get('ClinVarResult-Set').get('VariationArchive').get('@VariationID')


def clinvar_clinical_significance_pipeline(vus_df: pd.DataFrame) -> InternalResponse:
    # make a copy of the dataframe to be able to iterate through it whilst modifying the original dataframe
    new_vus_df = vus_df.copy()
    clinvar_variation_ids_dict = {}

    # get clinvar ids based on the rsid values
    for index, row in new_vus_df.iterrows():
        if row['RSID'] == "NORSID":
            current_app.logger.error(
                f"ClinVar clinical significance pipeline cannot run for variant without RSID!")
            vus_df.at[index, 'Clinvar error msg'] = "No RSID"
        else:
            current_app.logger.info(
                f"Retrieving information for:\n\tGene: {row['Gene']}\n\tChromosome: {row['Chr']}\n\tChromosome position: "
                f"{row['Position']}\n\tGenotype: {row['Genotype']}")

            retrieve_clinvar_ids_res = retrieve_clinvar_ids(row['RSID'])

            if retrieve_clinvar_ids_res.status != 200:
                vus_df.at[index, 'Clinvar error msg'] = "Failed to retrieve Clinvar id"
                current_app.logger.error(f"Retrieval of ClinVar id for {row['RSID']} failed 500!")
            else:
                var_clinvar_ids = retrieve_clinvar_ids_res.data

                # check the number of ClinVar IDs returned for a variant search
                if len(var_clinvar_ids) == 0:
                    vus_df.at[index, 'Clinvar error msg'] = "ClinVar ID has not been found!"
                elif len(var_clinvar_ids) > 1:
                    vus_df.at[index, 'Clinvar error msg'] = f"The following ClinVar IDs returned: {var_clinvar_ids}"
                # if only a single ClinVar ID has been returned
                else:
                    clinvar_variation_ids_dict[index] = var_clinvar_ids[0]

    # get clinvar information
    clinvar_ids = list(clinvar_variation_ids_dict.values())
    retrieve_clinvar_dict_res = retrieve_multiple_clinvar_dict(clinvar_ids)

    if retrieve_clinvar_dict_res.status != 200:
        current_app.logger.error(
            f"Retrieval of ClinVar document summaries failed 500!")
        return InternalResponse(None, 500)
    else:
        clinvar_dict = {}

        if retrieve_clinvar_dict_res.data is not None and retrieve_clinvar_dict_res.data.get('ClinVarResult-Set') is not None:
            variation_archive = retrieve_clinvar_dict_res.data.get('ClinVarResult-Set').get('VariationArchive')

            for v in variation_archive:
                clinvar_dict[v.get('@VariationID')] = v

            for index in clinvar_variation_ids_dict.keys():
                clinvar_variation_id = clinvar_variation_ids_dict[index]
                clinvar_variant_dict = clinvar_dict[clinvar_variation_id]

                classified_record_dict = clinvar_variant_dict.get('ClassifiedRecord')

                are_equivalent, error_msg = compare_clinvar_variant_with_expected_variant('GRCh37',
                                                                                          classified_record_dict,
                                                                                          vus_df.at[index, 'Gene'],  vus_df.at[index, 'Chr'],
                                                                                          vus_df.at[index, 'Position'])

                clinical_significance = extract_clinvar_germline_classification(classified_record_dict)

                if clinical_significance:
                    vus_df.at[index, 'Clinvar classification'] = clinical_significance['description']
                    vus_df.at[index, 'Clinvar classification last eval'] = clinical_significance['last_evaluated']
                    vus_df.at[index, 'Clinvar classification review status'] = clinical_significance['review_status']

                vus_df.at[index, 'Clinvar error msg'] = error_msg
                vus_df.at[index, 'Clinvar canonical spdi'] = extract_clinvar_canonical_spdi(classified_record_dict)
                vus_df.at[index, 'Clinvar variation id'] = clinvar_variation_id

    return InternalResponse(vus_df, 200)


def clinvar_clinical_significance_pipeline_single(genome_version: str, rsid: str, gene: str, chr: str,
                                            chr_pos: str) -> InternalResponse:
    is_success = True
    clinical_significance = {}
    canonical_spdi = ''
    variation_id = ''
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
                                                                                              clinvar_dict.get('ClinVarResult-Set').get('VariationArchive').get('ClassifiedRecord'),
                                                                                              gene, chr, chr_pos)
                    # if are_equivalent:
                    clinical_significance = extract_clinvar_germline_classification(clinvar_dict.get('ClinVarResult-Set').get('VariationArchive').get('ClassifiedRecord'))
                    canonical_spdi = extract_clinvar_canonical_spdi(clinvar_dict.get('ClinVarResult-Set').get('VariationArchive').get('ClassifiedRecord'))
                    variation_id = extract_clinvar_variation_id(clinvar_dict)
                    # else:
                    #     is_success = False

        return InternalResponse((is_success, clinical_significance, canonical_spdi, variation_id, error_msg), 200)


# TODO: retrieve other vital infor from clinvar (such as last modfied)
# Retrieve Clinvar variant classifications for every variant attempt to retrieve a corresponding
# ClinVar variant and extract its clinical significance.
def retrieve_clinvar_variant_classifications(vus_df: pd.DataFrame) -> InternalResponse:
    # TODO: add fix for when multiple rsids are found
    clinvar_clinical_significance_pipeline_res = clinvar_clinical_significance_pipeline(vus_df)

    if clinvar_clinical_significance_pipeline_res.status != 200:
        current_app.logger.error(
            f"ClinVar clinical significance pipeline failed!")
        return InternalResponse(None, 500)
    else:
        vus_df = clinvar_clinical_significance_pipeline_res.data

    return InternalResponse(vus_df, 200)


def get_last_saved_clinvar_update(clinvar_id: int) -> (int, str, str, str):
    clinvar_eval_date: AutoClinvarEvalDates = db.session.query(AutoClinvarEvalDates).filter(
        AutoClinvarEvalDates.clinvar_id == clinvar_id,
        AutoClinvarEvalDates.auto_clinvar_update_id.is_not(None)).order_by(
        desc(AutoClinvarEvalDates.eval_date)).first()

    auto_clinvar_update = clinvar_eval_date.auto_clinvar_update

    if auto_clinvar_update.last_evaluated is not None:
        clinvar_last_evaluated = datetime.strftime(auto_clinvar_update.last_evaluated, '%Y/%m/%d %H:%M')
    else:
        clinvar_last_evaluated = None

    return auto_clinvar_update.id, auto_clinvar_update.review_status, auto_clinvar_update.classification, clinvar_last_evaluated


def store_clinvar_info(clinvar_id: int, classification: str | None, review_status: str | None, last_eval: str | None,
                       is_new_vus: bool) -> \
        Tuple[bool, str]:
    clinvar_last_evaluated = None
    last_saved_classification = None

    if last_eval is not None and len(last_eval):
        clinvar_last_evaluated = datetime.strptime(last_eval, '%Y/%m/%d %H:%M')

    create_new_clinvar_update = True
    new_clinvar_update_id = None

    if not is_new_vus:
        auto_clinvar_update_id, last_saved_review_status, last_saved_classification, last_saved_eval_date = get_last_saved_clinvar_update(
            clinvar_id)

        # compare it to current clinvar info
        if last_eval == last_saved_eval_date and classification == last_saved_classification and review_status == last_saved_review_status:
            create_new_clinvar_update = False

    # create new clinvar update
    if is_new_vus or create_new_clinvar_update:
        new_clinvar_update = AutoClinvarUpdates(classification=classification,
                                                review_status=review_status,
                                                last_evaluated=clinvar_last_evaluated)
        db.session.add(new_clinvar_update)
        db.session.flush()

        new_clinvar_update_id = new_clinvar_update.id

    new_clinvar_eval_date = AutoClinvarEvalDates(eval_date=datetime.now(), clinvar_id=clinvar_id,
                                                 auto_clinvar_update_id=new_clinvar_update_id)
    db.session.add(new_clinvar_eval_date)

    return create_new_clinvar_update, last_saved_classification


# checks for the latest Clinvar updates
def get_updated_external_references_for_existing_vus(existing_vus_df: pd.DataFrame) -> InternalResponse:
    vus_df_copy = existing_vus_df.copy()

    # store the variants whose clinvar was updated
    updated_clinvar = []

    # load Clinvar variation ids
    for index, row in vus_df_copy.iterrows():
        is_variant_upload = 'Variant Id' in row.keys()

        if is_variant_upload:
            variant_id = row['Variant Id']
        else:
            variant_id = row['id']

        external_ref: ExternalReferences = db.session.query(ExternalReferences).filter(
            ExternalReferences.variant_id == variant_id, ExternalReferences.db_type == 'clinvar').one_or_none()

        # if variant has clinvar entry
        if external_ref is not None:
            clinvar = db.session.query(Clinvar).filter(Clinvar.external_clinvar_id == external_ref.id).one()

            retrieve_clinvar_dict_res = retrieve_clinvar_dict(clinvar.variation_id)

            if retrieve_clinvar_dict_res.status != 200:
                current_app.logger.error(
                    f"Retrieval of ClinVar document summary for ClinVar Variation Id {clinvar.variation_id} failed 500!")
                return InternalResponse(None, 500)
            else:
                clinvar_dict = retrieve_clinvar_dict_res.data

                # get latest clinvar classification
                clinvar_germline_classification = extract_clinvar_germline_classification(clinvar_dict.get('ClinVarResult-Set').get('VariationArchive').get('ClassifiedRecord'))
                latest_germline_classification = clinvar_germline_classification.get('description')

                is_clinvar_update_created, last_saved_classification = store_clinvar_info(clinvar.id,
                                                                                          latest_germline_classification,
                                                                                          clinvar_germline_classification.get(
                                                                                              'review_status'),
                                                                                          clinvar_germline_classification.get(
                                                                                              'last_evaluated'),
                                                                                          False)

                if is_clinvar_update_created:
                    if is_variant_upload:
                        chr = row['Chr']
                        pos = row['Position']
                        gene = row['Gene']
                        alt = row['Alt']
                        ref = row['Reference']
                    else:
                        chr = row['chromosome']
                        pos = row['chromosomePosition']
                        gene = row['gene']
                        alt = row['altAllele']
                        ref = row['refAllele']

                    updated_clinvar.append({'variant_id': variant_id, 'chromosome': chr,
                                            'chromosome_position': pos,
                                            'gene': gene, 'alt_allele': alt,
                                            'ref_allele': ref,
                                            'prev_classification': last_saved_classification,
                                            'new_classification': latest_germline_classification})

                # update Clinvar classification
                existing_vus_df.at[index, 'Clinvar classification'] = latest_germline_classification

    return InternalResponse({'existing_vus_df': existing_vus_df, 'clinvar_updates': updated_clinvar}, 200)


def scheduled_clinvar_updates():
    variants: List[Variants] = db.session.query(Variants).all()

    variant_summary = [get_variant_summary(v) for v in variants]

    vus_df = pd.DataFrame(variant_summary)

    get_updated_external_references_for_existing_vus_res = get_updated_external_references_for_existing_vus(vus_df)

    if get_updated_external_references_for_existing_vus_res.status != 200:
        current_app.logger.error(f'Failed to get Clinvar info')
        return InternalResponse({'isSuccess': False}, 500)
    else:
        clinvar_updates = get_updated_external_references_for_existing_vus_res.data['clinvar_updates']

        if len(clinvar_updates) > 0:
            # get all scientific members of staff
            recipients: List[ScientificMembers] = db.session.query(ScientificMembers).all()
            recipient_emails = [r.email for r in recipients]

            mail_message = Message(
                subject='Clinvar Classification Updates',
                bcc=recipient_emails,
                html=render_template('clinvar_update_email_template.html', updates=clinvar_updates,
                                     domain=os.environ.get('FRONT_URL')),
            )

            try:
                current_app.logger.info('Sending email with Clinvar updates.')
                mail.send(mail_message)
            except Exception as e:
                current_app.logger.error(f'Clinvar updates email was not sent successfully: {e}')

        try:
            # Commit the session to persist changes to the database
            db.session.commit()
            return InternalResponse({'isSuccess': True}, 200)
        except SQLAlchemyError as e:
            # Changes were rolled back due to an error
            db.session.rollback()

            current_app.logger.error(
                f'Rollback carried out since insertion of VariantsAcmgRules entry in DB failed due to error: {e}')
            return InternalResponse({'isSuccess': False}, 500)


def get_variant_clinvar_updates(clinvar_id: str):
    clinvar_updates_list = []

    clinvar: Clinvar = db.session.get(Clinvar, int(clinvar_id))

    eval_dates: List[AutoClinvarEvalDates] = clinvar.auto_clinvar_eval_dates

    # reversed to get dates in desc order
    eval_dates.reverse()

    for eval_date in eval_dates:
        update = None
        if eval_date.auto_clinvar_update_id is not None:
            auto_clinvar_update: AutoClinvarUpdates = eval_date.auto_clinvar_update
            update = {'classification': auto_clinvar_update.classification,
                      'reviewStatus': auto_clinvar_update.review_status,
                      'lastEval': datetime.strftime(auto_clinvar_update.last_evaluated, '%Y/%m/%d %H:%M')}

        clinvar_updates_list.append(
            {'dateChecked': datetime.strftime(eval_date.eval_date, '%d/%m/%Y %H:%M'), 'update': update})

    dates_with_updates = list(
        set([u['dateChecked'].split(" ")[0] for u in clinvar_updates_list if u['update'] is not None]))

    return Response(
        json.dumps({'isSuccess': True, 'clinvarUpdates': clinvar_updates_list, 'datesWithUpdates': dates_with_updates}),
        200, mimetype='application/json')


def get_last_clinvar_auto_update_date() -> str | None:
    auto_clinvar_eval_date: AutoClinvarEvalDates = db.session.query(AutoClinvarEvalDates).order_by(desc(AutoClinvarEvalDates.eval_date)).first()

    if auto_clinvar_eval_date is not None:
        return datetime.strftime(auto_clinvar_eval_date.eval_date, '%d/%m/%Y %H:%M')
    else:
        return None
