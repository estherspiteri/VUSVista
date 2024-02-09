import urllib
from urllib.parse import urlencode

import requests
from flask import Blueprint, Response, current_app, request
import json
import pandas as pd
from server.responses.internal_response import InternalResponse
from server.services.view_vus_service import retrieve_all_vus_from_db
from server.services.vus_preprocess_service import handle_vus_file

vus_views = Blueprint('vus_views', __name__)


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

    sample_phenotype_selection = request.form['samplePhenotypeSelection']

    # Parse the JSON string into a Python object
    if sample_phenotype_selection:
        sample_phenotype_selection_object = json.loads(sample_phenotype_selection)
    else:
        sample_phenotype_selection_object = []

    return handle_vus_file(file, sample_phenotype_selection_object, multiple_genes_selection_object)


@vus_views.route('/view', methods=['GET'])
def view_all_vus():
    current_app.logger.info(f"User requested to view all VUS")

    var_list = retrieve_all_vus_from_db()

    return Response(json.dumps({'isSuccess': True, 'vusList': var_list}), 200, mimetype='application/json')


@vus_views.route('/phenotype/<string:phenotype>', methods=['GET'])
def get_phenotype_terms(phenotype: str):
    url_encoded_phenotype = urllib.parse.quote(phenotype)
    url = f"https://hpo.jax.org/api/hpo/search?q={url_encoded_phenotype}&max=30&category=terms"

    hpo_res = requests.get(url)

    if hpo_res.status_code != 200:
        current_app.logger.error(
            f'Response failure {hpo_res.status_code}: {hpo_res.reason}')
        return Response(json.dumps({'isSuccess': False, 'hpoTerms': None}), 500, mimetype='application/json')
    else:
        terms = [{'ontologyId': x['ontologyId'],  'name': x['name']} for x in hpo_res.json()['terms']]
        return Response(json.dumps({'isSuccess': True, 'hpoTerms': terms}), 200, mimetype='application/json')