import json

from flask import Blueprint, current_app, Response, request

from server.services.phenotype_service import add_phenotype_to_existing_sample, remove_phenotype_to_sample, get_hpo_terms
from server.services.samples_service import retrieve_sample_from_db, retrieve_all_samples_from_db, delete_sample_entry, \
    update_variant_sample_hgvs

sample_views = Blueprint('sample_views', __name__)


@sample_views.route('/view/<string:sample_id>', methods=['GET'])
def get_sample(sample_id: str):
    current_app.logger.info(f"User requested to view sample with id: {sample_id}")

    res = retrieve_sample_from_db(sample_id)

    return Response(json.dumps({'isSuccess': res.data['isSuccess'], 'sample': res.data['sample_dict']}), 200,
                    mimetype='application/json')


@sample_views.route('/view', methods=['GET'])
def view_all_samples():
    current_app.logger.info(f"User requested to view all samples")

    sample_list = retrieve_all_samples_from_db()

    return Response(json.dumps({'isSuccess': True, 'sampleList': sample_list}), 200,
                    mimetype='application/json')


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


@sample_views.route('/delete/<string:sample_id>', methods=['DELETE'])
def delete_sample(sample_id: str):
    current_app.logger.info(f"Deleting sample with Id: {sample_id}")

    res = delete_sample_entry(sample_id)

    return Response(json.dumps({'isSuccess': res.status == 200}), res.status)


@sample_views.route('/edit-hgvs/<string:sample_id>/<string:variant_id>/<string:hgvs>', methods=['POST'])
def edit_variant_sample_hgvs(sample_id: str, variant_id: str, hgvs: str):
    current_app.logger.info(f"Editing variant-sample hgvs where variant Id: {variant_id}, sample Id: {sample_id} and new hgvs: {hgvs}")

    res = update_variant_sample_hgvs(sample_id, variant_id, hgvs)

    return Response(json.dumps({'isSuccess': res.status == 200}), res.status)
