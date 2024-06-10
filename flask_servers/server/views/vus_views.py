from datetime import datetime
from typing import List

import pandas as pd
from flask import Blueprint, Response, current_app, request
import json

from sqlalchemy.exc import SQLAlchemyError

from server import db
from server.models import GeneAttributes
from server.services.acmg_service import get_acmg_rules
from server.services.clinvar_service import get_variant_clinvar_updates
from server.services.view_vus_service import retrieve_all_vus_summaries_from_db, \
    retrieve_vus_from_db, delete_variant_entry
from server.services.vus_preprocess_service import handle_vus_file, handle_vus_from_form
from server.services.publications_service import add_publications_to_variant, \
    get_variant_publication_updates

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

    return get_variant_clinvar_updates(clinvar_id)


@vus_views.route('/get_publication_updates/<string:variant_id>', methods=['GET'])
def get_publication_updates(variant_id: str):
    current_app.logger.info(f"User requested all publication updates for Variant id {variant_id} ")

    return get_variant_publication_updates(variant_id)


@vus_views.route('/add_publications/<string:variant_id>', methods=['POST'])
def add_publications(variant_id: str):
    publication_links = request.form['publicationUrls']

    # Parse the JSON string into a Python object
    if publication_links:
        publication_links_list = json.loads(publication_links)
    else:
        publication_links_list = []

    return add_publications_to_variant(variant_id, publication_links_list)


@vus_views.route('/delete/<string:variant_id>', methods=['DELETE'])
def delete_variant(variant_id: str):
    current_app.logger.info(f"Deleting variant with Id: {variant_id}")

    res = delete_variant_entry(variant_id)

    return Response(json.dumps({'isSuccess': res.status == 200}), res.status)


