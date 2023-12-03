from flask import Blueprint, Response, current_app, request
import pandas as pd
import json

from server.helpers.data_helper import convert_df_to_list
from server.services.vus_preprocess_service import preprocess_vus, handle_vus_file

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

    # read the variant list into a dataframe
    # TODO: load from db
    try:
        var_df = pd.read_excel('final_vus.xlsx', header=0)
        var_df = var_df.fillna('')
    except Exception as e:
        current_app.logger.error(f'Error loading VUS from file: {e}')
        return Response(json.dumps({'isSuccess': False}), 500)

    # to match React camelCase syntax
    new_var_df = pd.DataFrame()
    new_var_df['vusID'] = var_df['VUS Id']
    new_var_df['chromosome'] = var_df['Chr']
    new_var_df['chromosomePosition'] = var_df['Position']
    new_var_df['gene'] = var_df['Gene']
    new_var_df['type'] = var_df['Type']
    new_var_df['genotype'] = var_df['Genotype']
    new_var_df['refAllele'] = var_df['Reference']
    new_var_df['observedAllele'] = var_df['Observed Allele']
    new_var_df['classification'] = var_df['Classification']
    new_var_df['rsid'] = var_df['RSID']
    new_var_df['rsidDbsnpVerified'] = var_df['RSID dbSNP verified']
    new_var_df['rsidDbsnpErrorMsgs'] = var_df['RSID dbSNP errorMsgs']
    new_var_df['clinvarErrorMsg'] = var_df['Clinvar error msg']
    new_var_df['clinvarClassification'] = var_df['Clinvar classification']
    new_var_df['clinvarClassificationLastEval'] = var_df['Clinvar classification last eval']
    new_var_df['clinvarClassificationReviewStatus'] = var_df['Clinvar classification review status']
    new_var_df['clinvarCanonicalSpdi'] = var_df['Clinvar canonical spdi']

    var_list = convert_df_to_list(new_var_df)

    return Response(json.dumps({'isSuccess': True, 'vusList': var_list}), 200, mimetype='application/json')

