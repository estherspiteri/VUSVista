from flask import Blueprint, Response, current_app, request
import pandas as pd
import json

from server.helpers.data_helper import convert_df_to_list, prep_vus_df_for_react
from server.services.vus_preprocess_service import handle_vus_file

vus_views = Blueprint('vus_views', __name__)


# @vus_views.route('/file/resume/<string:stage>', methods=['POST'])
# def resume_vus_preprocessing(stage: str):
#     current_app.logger.info(f"Resuming VUS file preprocessing at stage: {stage}")
#
#     return resume_vus_processing(stage)


# To include if you want VUS preprocessing, retrieval of RSIDs and retrieval of CLinVar classifications to be executed
# one after the other
@vus_views.route('/file', methods=['POST'])
def store_and_verify_vus_file():
    current_app.logger.info(f"User storing new VUS file")

    file = request.files['file']
    current_app.logger.info(f'Received file {file.filename} of type {file.content_type}')

    return handle_vus_file(file)


@vus_views.route('/view', methods=['GET'])
def view_all_vus():
    current_app.logger.info(f"User requested to view all VUS")

    # read the variant list into a dataframe
    # TODO: load from db
    try:
        var_df = pd.read_excel('final_vus.xlsx', header=0)
        var_df = var_df.fillna('')
    except Exception as e:
        current_app.logger.error(f'Error loading VUS from file: {e}')
        return Response(json.dumps({'isSuccess': False}), 500)

    new_var_df = prep_vus_df_for_react(var_df)

    var_list = convert_df_to_list(new_var_df)

    return Response(json.dumps({'isSuccess': True, 'vusList': var_list}), 200, mimetype='application/json')

