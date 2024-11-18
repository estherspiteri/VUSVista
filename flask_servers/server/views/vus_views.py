import threading
from datetime import datetime
from io import BytesIO
from typing import List

import pandas as pd
from flask import Blueprint, Response, current_app, request
import json

from flask_login import current_user
from sqlalchemy.exc import SQLAlchemyError

from server import db, file_upload_tasks
from server.models import GeneAttributes, FileUploadEvents, Variants
from server.services.acmg_service import get_acmg_rules
from server.services.clinvar_service import get_variant_clinvar_updates
from server.services.view_vus_service import retrieve_vus_summaries_from_db, \
    retrieve_vus_from_db, delete_variant_entry, add_samples_to_variant, add_new_sample_to_variant, \
    remove_sample_from_variant, update_variant_rsid, get_latest_added_vus
from server.services.vus_preprocess_service import handle_vus_file, handle_vus_from_form, check_for_multiple_genes, \
    check_for_existing_genes
from server.services.publications_service import add_publications_to_variant, \
    get_variant_publication_updates, remove_publications_to_variant

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


@vus_views.route('/file/multiple-genes-check', methods=['POST'])
def check_file_multiple_genes():
    current_app.logger.info(f"User checking for multiple genes in new VUS file")

    file = request.files['file']

    vus_df = pd.read_excel(file, header=0)

    multiple_genes = check_for_multiple_genes(vus_df)

    return Response(json.dumps({'isSuccess': True, 'multipleGenes': multiple_genes}), 200, mimetype='application/json')


@vus_views.route('/file/existing-genes-check', methods=['POST'])
def check_file_existing_genes():
    current_app.logger.info(f"User checking that genes in new VUS file can be found in DB")

    file = request.files['file']

    multiple_genes_selection = request.form['multipleGenesSelection']

    # Parse the JSON string into a Python object
    if multiple_genes_selection:
        multiple_genes_selection_arr = json.loads(multiple_genes_selection)
    else:
        multiple_genes_selection_arr = []

    # handle multiple genes selection
    vus_df = pd.read_excel(file, header=0)  # TODO: check re header
    current_app.logger.info(f'Number of VUS found in file: {len(vus_df)}')

    if len(multiple_genes_selection_arr) > 0:
        for selection in multiple_genes_selection_arr:
            vus_df.at[int(selection['index']), 'Gene'] = selection['gene']

    genes_not_found_in_db = check_for_existing_genes(vus_df)

    return Response(json.dumps({'isSuccess': True, 'genesNotInDb': genes_not_found_in_db}), 200, mimetype='application/json')


@vus_views.route('/file', methods=['POST'])
def store_and_verify_vus_file():
    current_app.logger.info(f"User storing new VUS file")

    file = request.files['file']
    current_app.logger.info(f'Received file {file.filename} of type {file.content_type}')

    vus_df = pd.read_excel(file, header=0)  # TODO: check re header
    current_app.logger.info(f'Number of VUS found in file: {len(vus_df)}')

    # TODO: merge multiple genes selection with genes not found selection
    #  further improvement: store updated vus file in database after each stage

    # handle multiple genes selection
    multiple_genes_selection = request.form['multipleGenesSelection']

    # Parse the JSON string into a Python object
    if multiple_genes_selection:
        multiple_genes_selection_arr = json.loads(multiple_genes_selection)
    else:
        multiple_genes_selection_arr = []

    if len(multiple_genes_selection_arr) > 0:
        for selection in multiple_genes_selection_arr:
            vus_df.at[int(selection['index']), 'Gene'] = selection['gene']

    # handle genes not found selection
    genes_not_found_selection = request.form['genesNotFoundSelection']

    # Parse the JSON string into a Python object
    if genes_not_found_selection:
        genes_not_found_selection_arr = json.loads(genes_not_found_selection)
    else:
        genes_not_found_selection_arr = []

    if len(genes_not_found_selection_arr) > 0:
        for selection in genes_not_found_selection_arr:
            vus_df.at[int(selection['index']), 'Gene'] = selection['gene']

    # Create a BytesIO object
    output = BytesIO()

    # Write the DataFrame to the BytesIO object as an Excel file
    vus_df.to_excel(output, index=False, engine='openpyxl')

    # Get the binary content
    file_data = output.getvalue()

    file_upload_event: FileUploadEvents = FileUploadEvents(file_name=file.filename, file_data=file_data, date_created=datetime.now(), date_processed=None, scientific_members_id=current_user.id)
    db.session.add(file_upload_event)

    try:
        # Commit the session to persist changes to the database
        db.session.commit()

        task_id = file_upload_event.id

        file_upload_tasks[task_id] = {'isSuccess': None, 'taskId': task_id}

        # handle_vus_file(task_id, file, multiple_genes_selection_object)

        return Response(json.dumps({'isSuccess': True, 'taskId': task_id}), 200, mimetype='application/json')
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since storing file {file.filename}\'s upload event in DB failed due to error: {e}')
        return Response(json.dumps({'isSuccess': False, 'taskId': None}), 200, mimetype='application/json')


@vus_views.route('file/check-status/<string:task_ids>', methods=['GET'])
def check_file_upload_status(task_ids: str):
    current_app.logger.info(f"Checking upload status of files with task ids {task_ids}")

    task_ids_arr = [int(task_id) for task_id in task_ids.split(',')]

    statuses = []

    for task_id in task_ids_arr:
        statuses.append(file_upload_tasks.get(task_id, {'taskId': task_id, 'isSuccess': None}))

    return Response(json.dumps({'statuses': statuses}), 200, mimetype='application/json')


@vus_views.route('/view', methods=['GET'])
def view_all_vus():
    current_app.logger.info(f"User requested to view all VUS")

    variants: List[Variants] = db.session.query(Variants).all()
    var_list = retrieve_vus_summaries_from_db(variants)

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


@vus_views.route('/remove_publication/<string:variant_id>/<int:publication_id>', methods=['POST'])
def remove_publications(variant_id: str, publication_id: int):
    current_app.logger.info(f"Removing publication {publication_id} from variant with Id: {variant_id}")

    return remove_publications_to_variant(variant_id, publication_id)


@vus_views.route('/delete/<string:variant_id>', methods=['DELETE'])
def delete_variant(variant_id: str):
    current_app.logger.info(f"Deleting variant with Id: {variant_id}")

    res = delete_variant_entry(variant_id)

    return Response(json.dumps({'isSuccess': res.status == 200}), res.status)


@vus_views.route('/add-existing-samples/<string:variant_id>', methods=['POST'])
def add_variant_existing_samples(variant_id: str):
    samples_to_add = request.form['samplesToAdd']

    # Parse the JSON string into a Python object
    if samples_to_add:
        samples_to_add_list = json.loads(samples_to_add)
    else:
        samples_to_add_list = []

    current_app.logger.info(f"Add samples {[s['sampleId'] for s in samples_to_add_list]} to variant {variant_id}")

    res = add_samples_to_variant(int(variant_id), samples_to_add_list)

    return Response(json.dumps({'isSuccess': res.status == 200, 'updatedSamples': res.data['updatedSamples'], "updatedNotVariantSamples": res.data['updatedNotVariantSamples'], "updatedPhenotypes": res.data["updatedPhenotypes"]}), res.status)


@vus_views.route('/add-new-sample/<string:variant_id>', methods=['POST'])
def add_variant_new_sample(variant_id: str):
    sample_to_add = request.form['sampleToAdd']

    # Parse the JSON string into a Python object
    if sample_to_add:
        sample_to_add_dict = json.loads(sample_to_add)
    else:
        sample_to_add_dict = {}

    current_app.logger.info(f"Add new sample {sample_to_add_dict['sampleId']} to variant {variant_id}")

    res = add_new_sample_to_variant(int(variant_id), sample_to_add_dict)

    return Response(json.dumps({'isSuccess': res.status == 200, 'updatedSamples': res.data['updatedSamples'], "updatedNotVariantSamples": res.data['updatedNotVariantSamples'], "updatedPhenotypes": res.data["updatedPhenotypes"]}), res.status)


@vus_views.route('/remove-samples/<string:variant_id>', methods=['POST'])
def remove_samples(variant_id: str):
    sample_ids_to_remove = request.form['sampleIdsToRemove']

    # Parse the JSON string into a Python object
    if sample_ids_to_remove:
        sample_ids_to_remove_list = json.loads(sample_ids_to_remove)
    else:
        sample_ids_to_remove_list = {}

    current_app.logger.info(f"Remove samples {sample_ids_to_remove_list} from variant {variant_id}")

    res = remove_sample_from_variant(int(variant_id), sample_ids_to_remove_list)

    return Response(json.dumps({'isSuccess': res.status == 200, 'updatedSamples': res.data['updatedSamples'], "updatedNotVariantSamples": res.data['updatedNotVariantSamples'], "updatedPhenotypes": res.data["updatedPhenotypes"]}), res.status)


@vus_views.route('/update-rsid/<int:variant_id>/<string:new_rsid>', methods=['POST'])
def update_rsid(variant_id: int, new_rsid: str):
    current_app.logger.info(f"Updating RSID for variant ID {variant_id} to: {new_rsid}")

    res = update_variant_rsid(variant_id, new_rsid)

    return Response(json.dumps({'isSuccess': res.status == 200, 'updatedExternalRefData': res.data['updated_external_ref_data']}), res.status)


@vus_views.route('/latest-uploaded-vus/<int:number_of_variants>', methods=['GET'])
def get_latest_uploaded_vus(number_of_variants: int):
    current_app.logger.info(f"Retrieving last {number_of_variants} uploaded variants")

    res = get_latest_added_vus(number_of_variants)

    return Response(json.dumps({'isSuccess': True, 'vusList': res.data['var_list']}), res.status)
