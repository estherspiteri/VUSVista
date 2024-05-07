import json

from flask import Blueprint, current_app, Response, request

from server import db
from server.helpers.data_helper import alchemy_encoder, get_variant_summary
from server.models import Variants
from server.services.litvar_service import get_publications
from server.services.vus_publications_service import get_publications_by_variant_id_from_db

publication_views = Blueprint('publication_views', __name__)


# TODO: check if can be removed
@publication_views.route('/getByRsid/<string:rsid>', methods=['GET'])
def get_publications_of_variant_by_rsid(rsid: str):
    # TODO: check that rsid starts with 'rs'
    current_app.logger.info(f"User requested retrieval of publications for variant with rsid {rsid}")

    get_publications_res = get_publications(None, rsid, '')

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

        return Response(json.dumps({'isSuccess': True,
                                    'publicationSearch': {'publications': publication_list,
                                                          "isLitvarIdFound": len(get_publications_res.data) > 0}}),
                        get_publications_res.status)


@publication_views.route('/getByVariantId/<string:variant_id>', methods=['GET'])
def get_publications_of_variant_by_variant_id(variant_id: str):
    current_app.logger.info(f"User requested retrieval of publications for variant with variant id {variant_id}")

    variant, publications = get_publications_by_variant_id_from_db(variant_id)

    variant_summary = get_variant_summary(variant)

    return Response(json.dumps({'isSuccess': True, 'variant': variant_summary, 'publications': publications}), 200)


@publication_views.route('/getWithOptionalText/<string:variant_id>/<string:rsid>/<string:optional_text>', methods=['GET'])
def get_publications_of_variant_by_rsid_with_optional_text(variant_id: str, rsid: str, optional_text: str):
    current_app.logger.info(f"User requested retrieval of publications for variant with Id {variant_id}, RSID {rsid} and optional text {optional_text}")

    get_publications_res = get_publications(None, rsid, optional_text)

    publication_list = []

    for publication in get_publications_res.data:
        encoded_publication = alchemy_encoder(publication)

        # changing keys of dictionary
        encoded_publication['date'] = encoded_publication.pop('date_published')

        publication_list.append(encoded_publication)

    variant = db.session.query(Variants).filter(Variants.id == variant_id).first()

    variant_summary = get_variant_summary(variant)

    return Response(json.dumps({'isSuccess': get_publications_res.status == 200, 'variant': variant_summary, 'publications': publication_list}), get_publications_res.status)

