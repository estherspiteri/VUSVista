import json

from flask import Blueprint, current_app, Response, jsonify
from sqlalchemy.orm import DeclarativeMeta

from server.helpers.data_helper import alchemy_encoder
from server.services.litvar_service import get_publications

litvar_views = Blueprint('litvar_views', __name__)


@litvar_views.route('/publications/<string:rsid>', methods=['GET'])
def get_publications_of_variant(rsid: str):
    # TODO: check that rsid starts with 'rs'
    current_app.logger.info(f"User requested retrieval of publications for variant with rsid {rsid}")

    get_publications_res = get_publications(rsid)

    if get_publications_res.status != 200:
        current_app.logger.error(f'LitVar publication retrieval query failed 500')
        return Response(json.dumps({'isSuccess': False}), 500)
    else:
        publication_list = []

        for publication in get_publications_res.data:
            encoded_publication = alchemy_encoder(publication)

            # changing keys of dictionary
            encoded_publication['date'] = encoded_publication.pop('date_published')

            publication_list.append(encoded_publication)

        print(publication_list)
        return Response(json.dumps({'isSuccess': True,
                                    'publicationSearch': {'publications': publication_list,
                                                          "isLitvarIdFound": len(get_publications_res.data) > 0}}),
                        get_publications_res.status)
