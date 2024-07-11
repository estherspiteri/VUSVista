import json

from flask import Blueprint, current_app, Response

from server.services.clinvar_service import get_last_clinvar_auto_update_date
from server.services.publications_service import get_last_pub_auto_update_date
from server.services.view_vus_service import get_latest_added_vus

homepage_views = Blueprint('homepage_views', __name__)


@homepage_views.route('/home/<int:number_of_variants>', methods=['GET'])
def get_latest_uploaded_vus(number_of_variants: int):
    current_app.logger.info(f"Retrieving last {number_of_variants} uploaded variants")

    get_latest_added_vus_res = get_latest_added_vus(number_of_variants)

    last_clinvar_auto_update_date = get_last_clinvar_auto_update_date()
    last_pub_auto_update_date = get_last_pub_auto_update_date()

    return Response(json.dumps({'isSuccess': True, 'data': {'vusList': get_latest_added_vus_res.data['var_list'], 'lastClinvarUpdateDate': last_clinvar_auto_update_date, 'lastPubUpdateDate': last_pub_auto_update_date}}), 200)
