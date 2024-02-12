import json

from flask import Blueprint, current_app, Response

from server.services.litvar_service import get_publications

litvar_views = Blueprint('litvar_views', __name__)


@litvar_views.route('/publications/<string:rsid>', methods=['GET'])
def get_publications_of_variant(rsid: str):
    # TODO: check that rsid starts with 'rs'
    current_app.logger.info(f"User requested retrieval of publications for variant with rsid {rsid}")

    get_publications_res = get_publications(rsid)
    return Response(json.dumps(get_publications_res.data), get_publications_res.status)