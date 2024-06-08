from datetime import datetime
from typing import List

import pandas as pd
from flask import Blueprint, Response, current_app, request
import json

from sqlalchemy.exc import SQLAlchemyError

from server import db
from server.models import GeneAttributes, Clinvar, AutoClinvarEvalDates, AutoClinvarUpdates, Variants, \
    VariantsPublications, Publications, AutoPublicationEvalDates
from server.services.acmg_service import get_acmg_rules, add_acmg_rule_to_variant, remove_acmg_rule_from_variant
from server.services.view_vus_service import retrieve_all_vus_summaries_from_db, \
    retrieve_vus_from_db
from server.services.vus_preprocess_service import handle_vus_file, preprocess_vus, handle_vus_from_form
from server.services.vus_publications_service import get_publication_info, merge_2_sets_of_publications, \
    store_variant_publications_in_db, get_publications_by_variant_id_from_db

vus_views = Blueprint('vus_views', __name__)


@vus_views.route('/upload', methods=['POST'])
def store_vus():
    current_app.logger.info(f"User storing new VUS")

    vus = request.form['vus']

    # Parse the JSON string into a Python object
    if vus:
        vus_object = json.loads(vus)
        vus_object['samples'] = ','.join([str(elem) for elem in vus_object['samples']])
    else:
        vus_object = {}

    # Convert dictionary values into arrays containing the string value
    vus_dict = {key: [value] for key, value in vus_object.items()}

    vus_df = pd.DataFrame.from_dict(vus_dict)

    return handle_vus_from_form(vus_df)


@vus_views.route('/file', methods=['POST'])
def store_and_verify_vus_file():
    current_app.logger.info(f"User storing new VUS file")

    file = request.files['file']
    current_app.logger.info(f'Received file {file.filename} of type {file.content_type}')

    multiple_genes_selection = request.form['multipleGenesSelection']

    # Parse the JSON string into a Python object
    if multiple_genes_selection:
        multiple_genes_selection_object = json.loads(multiple_genes_selection)
    else:
        multiple_genes_selection_object = []

    return handle_vus_file(file, multiple_genes_selection_object)


@vus_views.route('/view', methods=['GET'])
def view_all_vus():
    current_app.logger.info(f"User requested to view all VUS")

    var_list = retrieve_all_vus_summaries_from_db()

    return Response(json.dumps({'isSuccess': True, 'vusList': var_list}), 200, mimetype='application/json')


@vus_views.route('/view/<string:vus_id>', methods=['GET'])
def get_vus(vus_id: int):
    current_app.logger.info(f"User requested to view vus with id {vus_id}")

    var_list = retrieve_vus_from_db(vus_id)

    acmg_rules = get_acmg_rules()

    return Response(json.dumps({'isSuccess': True, 'vus': var_list, 'acmgRules': acmg_rules}), 200,
                    mimetype='application/json')


@vus_views.route('/gene/<string:gene_name>', methods=['GET'])
def verify_gene(gene_name: str):
    current_app.logger.info(f"User requested to verify gene {gene_name}")

    # Retrieving gene that matches gene_name from Gene Attributes table
    gene_attribute: GeneAttributes = db.session.query(GeneAttributes).filter(GeneAttributes.attribute_name == 'gene_name', GeneAttributes.attribute_value == gene_name.upper()).one_or_none()

    gene_id = None
    if gene_attribute is not None:
        gene_id = gene_attribute.gene_id

    return Response(json.dumps({'isSuccess': gene_id is not None, 'geneId': gene_id}), 200, mimetype='application/json')


@vus_views.route('/all-acmg-rules', methods=['GET'])
def get_all_acmg_rules():
    current_app.logger.info(f"User requested all acmg rules ")

    acmg_rules = get_acmg_rules()

    return Response(json.dumps({'isSuccess': True, 'acmgRules': acmg_rules}), 200, mimetype='application/json')


@vus_views.route('/get_clinvar_updates/<string:clinvar_id>', methods=['GET'])
def get_clinvar_updates(clinvar_id: str):
    current_app.logger.info(f"User requested all clinvar updates for Clinvar id {clinvar_id} ")

    clinvar_updates_list = []

    clinvar: Clinvar = db.session.get(Clinvar, int(clinvar_id))

    eval_dates: List[AutoClinvarEvalDates] = clinvar.auto_clinvar_eval_dates

    # reversed to get dates in desc order
    eval_dates.reverse()

    for eval_date in eval_dates:
        update = None
        if eval_date.auto_clinvar_update_id is not None:
            auto_clinvar_update: AutoClinvarUpdates = eval_date.auto_clinvar_update
            update = {'classification': auto_clinvar_update.classification, 'reviewStatus': auto_clinvar_update.review_status, 'lastEval': datetime.strftime(auto_clinvar_update.last_evaluated, '%Y/%m/%d %H:%M')}

        clinvar_updates_list.append({'dateChecked': datetime.strftime(eval_date.eval_date, '%d/%m/%Y %H:%M'), 'update': update})

    dates_with_updates = list(set([u['dateChecked'].split(" ")[0] for u in clinvar_updates_list if u['update'] is not None]))

    return Response(json.dumps({'isSuccess': True, 'clinvarUpdates': clinvar_updates_list, 'datesWithUpdates': dates_with_updates}), 200, mimetype='application/json')


@vus_views.route('/get_publication_updates/<string:variant_id>', methods=['GET'])
def get_publication_updates(variant_id: str):
    current_app.logger.info(f"User requested all publication updates for Variant id {variant_id} ")

    variant: Variants = db.session.get(Variants, int(variant_id))

    variant_publication_entries: List[VariantsPublications] = variant.variants_publications

    auto_publication_updates = {}
    manual_publication_updates = {}
    for vp in variant_publication_entries:
        publication: Publications = db.session.query(Publications).get(vp.publication_id)
        date = datetime.strftime(vp.date_added, '%d/%m/%Y %H:%M')

        if vp.is_manually_added:
            update_arr = manual_publication_updates.get(date, [])
            manual_publication_updates[date] = update_arr + [{"title": publication.title, "doi": publication.doi, "link": publication.link, "isManuallyAdded": True}]
        else:
            update_arr = auto_publication_updates.get(date, [])
            auto_publication_updates[date] = update_arr + [{"title": publication.title, "doi": publication.doi, "link": publication.link, "isManuallyAdded": False}]

    auto_pub_eval_dates: List[AutoPublicationEvalDates] = db.session.query(AutoPublicationEvalDates).filter(AutoPublicationEvalDates.variant_id == variant_id).all()

    # get both manual and automatic updates dates
    all_auto_pub_eval_dates = [d.eval_date for d in auto_pub_eval_dates]
    all_manual_pub_eval_dates = [vp.date_added for vp in variant_publication_entries if vp.is_manually_added == True]
    all_dates = all_auto_pub_eval_dates + all_manual_pub_eval_dates

    # reversed to get dates in desc order
    all_dates.sort()
    all_dates.reverse()

    all_dates_str = [datetime.strftime(d, '%d/%m/%Y %H:%M') for d in all_dates]
    all_dates_str_unique = []
    for d in all_dates_str:
        if d not in all_dates_str_unique:
            all_dates_str_unique.append(d)

    all_auto_pub_eval_dates_str = [datetime.strftime(d.eval_date, '%d/%m/%Y %H:%M') for d in auto_pub_eval_dates]

    variant_publications_update_list = []
    for date in all_dates_str_unique:
        if date in all_auto_pub_eval_dates_str and date not in auto_publication_updates.keys():
            pub_updates = None

        else:
            pub_updates = []

            if date in auto_publication_updates.keys():
                pub_updates += auto_publication_updates[date]

            if date in manual_publication_updates.keys():
                pub_updates += manual_publication_updates[date]

        variant_publications_update_list.append({'lastEval': date,
                                                 "publicationUpdates": pub_updates})

    dates_with_updates = list(set([u['lastEval'].split(" ")[0] for u in variant_publications_update_list if u['publicationUpdates'] is not None]))

    return Response(json.dumps({'isSuccess': True, 'variantPublicationUpdates': variant_publications_update_list, 'datesWithUpdates': dates_with_updates}), 200, mimetype='application/json')


@vus_views.route('/add_publications/<string:variant_id>', methods=['POST'])
def add_publications(variant_id: str):
    publication_links = request.form['publicationUrls']

    # Parse the JSON string into a Python object
    if publication_links:
        publication_links_list = json.loads(publication_links)
    else:
        publication_links_list = []

    variant: Variants = db.session.query(Variants).get(variant_id)
    variant_publication_ids = [vp.publication_id for vp in variant.variants_publications]
    variant_publications: List[Publications] = db.session.query(Publications).filter(Publications.id.in_(variant_publication_ids)).all()

    added_publications: List[Publications] = []
    for link in publication_links_list:
        link_publication_res = get_publication_info(link)

        if link_publication_res.status != 200:
            current_app.logger.error(
                f'Retrieval of information for the user provided literature link failed 500')
        else:
            new_pub = link_publication_res.data
            added_publications.append(new_pub)

    final_pub_list = merge_2_sets_of_publications(added_publications, variant_publications)

    # extract the publications not yet included for the variant
    pub_not_in_variant = [p for p in final_pub_list if p.id not in variant_publication_ids]
    pub_not_in_variant_dois = [p.doi for p in pub_not_in_variant]
    store_variant_publications_in_db(pub_not_in_variant, variant.id, pub_not_in_variant_dois)

    try:
        # Commit the session to persist changes to the database
        db.session.commit()

        variant, publications = get_publications_by_variant_id_from_db(variant_id)
        return Response(json.dumps({'isSuccess': True, 'publications': publications}), 200, mimetype='application/json')
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since insertion of new publication entries in DB failed due to error: {e}')
        return Response(json.dumps({'isSuccess': False, "publications": None}), 200, mimetype='application/json')


