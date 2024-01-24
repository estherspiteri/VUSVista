import json

from flask import Blueprint, current_app, Response

from server.services.view_samples_service import retrieve_all_samples_from_db

sample_views = Blueprint('sample_views', __name__)


@sample_views.route('/view', methods=['GET'])
def view_all_samples():
    current_app.logger.info(f"User requested to view all samples")

    sample_list = retrieve_all_samples_from_db()

    return Response(json.dumps({'isSuccess': True, 'sampleList': sample_list}), 200,
                    mimetype='application/json')

