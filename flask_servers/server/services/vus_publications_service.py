import json
from datetime import datetime
from typing import List, Dict
from urllib.parse import urlencode

import pandas as pd
import requests
from flask import current_app

from server import db
from server.helpers.data_helper import alchemy_encoder
from server.helpers.db_access_helper import get_variant_from_db
from server.models import Publications, Variants
from server.responses.internal_response import InternalResponse
from server.services.litvar_service import get_litvar_publications, get_publications


def get_publication_info(publication_link: str) -> InternalResponse:
    # using Open Access Button: https://openaccessbutton.org/api
    url_encoded_link = urlencode({'id': publication_link})
    url = f"https://bg.api.oa.works/find?{url_encoded_link}"

    open_access_res = requests.get(url)

    if open_access_res.status_code != 200:
        current_app.logger.error(
            f'Response failure {open_access_res.status_code}: {open_access_res.reason}')
        return InternalResponse(None, open_access_res.status_code, open_access_res.reason)
    else:
        metadata = open_access_res.json()['metadata']

        date = None
        if 'published' in metadata.keys():
            date = metadata['published']
        elif 'year' in metadata.keys():
            date = metadata['year'] + '-01-01' # TODO: fix - misleading on front-end

        if date is not None:
            date = datetime.strptime(date, '%Y-%m-%d')

        # eliminate any keywords included
        updated_doi = metadata.get('doi', None)
        if updated_doi is not None:
            updated_doi = updated_doi.replace("\"", "").split(',')[0]

        publication = Publications(title=metadata.get('title', None), pmid=metadata.get('pmid', None),
                                   doi=updated_doi,
                                   abstract=metadata.get('abstract', None),
                                   date_published=date,
                                   link=publication_link)

        return InternalResponse(publication, 200)


def extract_publications_already_stored_in_db(pub_list: List[Publications]) -> (List[Publications], List[Publications]):
    doi_list = [p.doi for p in pub_list]

    matching_db_publications: List[Publications] = db.session.query(Publications).filter(Publications.doi.in_(doi_list)).all()
    matching_db_dois = [p.doi for p in matching_db_publications]

    new_pub_list = [p for p in pub_list if p.doi not in matching_db_dois]

    return new_pub_list, matching_db_publications


# merge user provided publications with those from LitVar for new variants
def merge_user_and_litvar_publications(user_pub_list: List[Publications], litvar_publications: List[Publications]) -> \
        (List[Publications], List[Publications]):
    litvar_publications_dois = [p.doi for p in litvar_publications]

    # check if a publication is listed in both the user provided publication links and LitVar's publications
    # by comparing the doi
    common_publication_dois = [p.doi for p in user_pub_list if p.doi in litvar_publications_dois]

    if len(common_publication_dois) == 0:
        final_pub_list = user_pub_list + litvar_publications
    else:
        # remove the publication provided by the user
        updated_user_pub = [p for p in user_pub_list if p.doi not in common_publication_dois]
        final_pub_list = updated_user_pub + litvar_publications

    final_pub_list, existing_pub_list = extract_publications_already_stored_in_db(final_pub_list)

    return final_pub_list, existing_pub_list


# merge user provided publications with those already stored in the database
def merge_user_and_db_publications(user_pub_list: List[Publications], db_publications: List[Publications]) -> \
        (List[Publications], List[Publications]):
    db_publications_dois = [p.doi for p in db_publications]

    # check if a publication is listed in both the user provided publication links and the publications already stored
    # in the db by comparing the doi
    common_publication_dois = [p.doi for p in user_pub_list if p.doi in db_publications_dois]

    if len(common_publication_dois) == 0:
        final_pub_list = user_pub_list
    else:
        # remove the publications already stored in the db
        updated_user_pub = [p for p in user_pub_list if p.doi not in common_publication_dois]
        final_pub_list = updated_user_pub

    final_pub_list, existing_pub_list = extract_publications_already_stored_in_db(final_pub_list)

    return final_pub_list, existing_pub_list


def retrieve_and_store_variant_publications(vus_df: pd.DataFrame, variants_already_stored_in_db: bool):
    current_app.logger.info('Retrieving publications for variants')

    # iterate through the dataframe
    for index, row in vus_df.iterrows():
        publication_links: List[Publications] = []
        litvar_publications: List[Publications] = []

        # check if user provided any literature links
        links = str(row['Literature Links']).replace(' ', '')

        if len(links) > 0:
            links_arr = links.split('|')
        else:
            links_arr = []

        for link in links_arr:
            link_publication_res = get_publication_info(link)

            if link_publication_res.status != 200:
                current_app.logger.error(
                    f'Retrieval of information for the user provided literature link failed 500')
                return InternalResponse(None, 500)
            else:
                publication_links.append(link_publication_res.data)

        if not variants_already_stored_in_db and len(row['RSID']) > 0 and row['RSID'] != 'NORSID':
            # retrieve LitVar publications
            litvar_publications_res = get_publications(row['RSID'])

            if litvar_publications_res.status != 200:
                current_app.logger.error(
                    f'Retrieval of information for the user provided literature link failed 500')
                return InternalResponse(None, 500)
            else:
                litvar_publications = litvar_publications_res.data

        # retrieve variant
        variant = get_variant_from_db(row)

        if variants_already_stored_in_db:
            vus_new_publications, vus_existing_publications = merge_user_and_db_publications(publication_links,
                                                                                  variant.publications)
        else:
            vus_new_publications, vus_existing_publications = merge_user_and_litvar_publications(publication_links,
                                                                                      litvar_publications)

        # store new publications
        for p in vus_new_publications:
            db.session.add(p)

            # set up relationship between the variant & its publications
            variant.publications.append(p)

        # set up relationships of variants with existing publications
        for p in vus_existing_publications:
            variant.publications.append(p)

    return InternalResponse(None, 200)


def get_publications_by_variant_id(variant_id: str) -> (Variants, List[Dict]):
    variant: Variants = db.session.query(Variants).filter(Variants.id == variant_id).one_or_none()

    publication_list = []

    var_publications: List[Publications] = variant.publications

    for p in var_publications:
        encoded_publication = alchemy_encoder(p)

        # changing keys of dictionary
        encoded_publication['publicationId'] = encoded_publication.pop('id')
        encoded_publication['date'] = encoded_publication.pop('date_published').strftime('%Y/%m/%d')
        encoded_publication['isSupplementaryMaterialMatch'] = encoded_publication.pop('match_in_sup_material')

        if p.authors is not None:
            encoded_publication['authors'] = p.authors.replace('\"', '').strip('{}').split(',')

        publication_list.append(encoded_publication)

    return variant, publication_list
