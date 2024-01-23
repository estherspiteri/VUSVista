from flask import Blueprint, Response, current_app, request
import json

from server.services.view_vus_service import retrieve_all_vus_from_db
from server.services.vus_preprocess_service import handle_vus_file

vus_views = Blueprint('vus_views', __name__)


@vus_views.route('/file', methods=['POST'])
def store_and_verify_vus_file():
    current_app.logger.info(f"User storing new VUS file")

    file = request.files['file']
    current_app.logger.info(f'Received file {file.filename} of type {file.content_type}')

    return handle_vus_file(file)


@vus_views.route('/view', methods=['GET'])
def view_all_vus():
    current_app.logger.info(f"User requested to view all VUS")

    var_list = retrieve_all_vus_from_db()

    return Response(json.dumps({'isSuccess': True, 'vusList': var_list}), 200, mimetype='application/json')

