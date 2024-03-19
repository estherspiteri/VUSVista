import json
import urllib
from urllib.parse import urlencode

import requests
from flask import Blueprint, current_app, Response, request

from server.services.acmg_service import get_acmg_rules, add_acmg_rule_to_sample_variant, \
    remove_acmg_rule_to_sample_variant
from server.services.phenotype_service import add_phenotype_to_existing_sample, remove_phenotype_to_sample, get_hpo_terms
from server.services.view_samples_service import retrieve_all_samples_from_db, retrieve_sample_from_db

sample_views = Blueprint('sample_views', __name__)


@sample_views.route('/view/<string:sample_id>', methods=['GET'])
def get_sample(sample_id: str):
    current_app.logger.info(f"User requested to view sample with id: {sample_id}")

    sample = retrieve_sample_from_db(sample_id)

    acmg_rules = get_acmg_rules()

    return Response(json.dumps({'isSuccess': True, 'sample': sample, 'acmgRules': acmg_rules}), 200,
                    mimetype='application/json')


@sample_views.route('/view', methods=['GET'])
def view_all_samples():
    current_app.logger.info(f"User requested to view all samples")

    sample_list = retrieve_all_samples_from_db()

    return Response(json.dumps({'isSuccess': True, 'sampleList': sample_list}), 200,
                    mimetype='application/json')


@sample_views.route('/add-acmg-rule', methods=['POST'])
def add_acmg_rule():
    current_app.logger.info(f"Adding ACMG rule")

    sample_id = request.form['sampleId']
    variant_id = request.form['variantId']
    rule_id = request.form['ruleId']

    res = add_acmg_rule_to_sample_variant(sample_id, int(variant_id), rule_id)

    return Response(json.dumps({'isSuccess': res.status == 200}), res.status)


@sample_views.route('/remove-acmg-rule', methods=['POST'])
def remove_acmg_rule():
    current_app.logger.info(f"Removing ACMG rule")

    sample_id = request.form['sampleId']
    variant_id = request.form['variantId']
    rule_id = request.form['ruleId']

    res = remove_acmg_rule_to_sample_variant(sample_id, int(variant_id), rule_id)

    return Response(json.dumps({'isSuccess': res.status == 200}), res.status)


@sample_views.route('/phenotype/<string:phenotype>', methods=['GET'])
def get_phenotype_terms(phenotype: str):
    hpo_res = get_hpo_terms(phenotype)

    return Response(json.dumps(hpo_res.data), hpo_res.status, mimetype='application/json')


@sample_views.route('/add-phenotype', methods=['POST'])
def add_phenotype():
    current_app.logger.info(f"Adding Phenotype")

    sample_id = request.form['sampleId']
    phenotype = request.form['phenotype']

    # Parse the JSON string into a Python object
    if phenotype:
        phenotype_dict = json.loads(phenotype)
    else:
        phenotype_dict = {}

    res = add_phenotype_to_existing_sample(sample_id, phenotype_dict)

    return Response(json.dumps({'isSuccess': res.status == 200}), res.status)


@sample_views.route('/remove-phenotype', methods=['POST'])
def remove_phenotype():
    current_app.logger.info(f"Removing Phenotype")

    sample_id = request.form['sampleId']
    phenotype = request.form['phenotype']

    # Parse the JSON string into a Python object
    if phenotype:
        phenotype_dict = json.loads(phenotype)
    else:
        phenotype_dict = {}

    res = remove_phenotype_to_sample(sample_id, phenotype_dict)

    return Response(json.dumps({'isSuccess': res.status == 200}), res.status)
