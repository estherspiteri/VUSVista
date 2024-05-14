import json
from datetime import datetime
from typing import List, Dict
from urllib.parse import urlencode

import pandas as pd
import requests
from flask import current_app
from sqlalchemy.exc import SQLAlchemyError

from server import db
from server.helpers.data_helper import alchemy_encoder
from server.helpers.db_access_helper import get_variant_from_db
from server.models import Publications, Variants, ExternalReferences
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
            date = metadata['year'] + '-01-01'  # TODO: fix - misleading on front-end

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

    matching_db_publications: List[Publications] = db.session.query(Publications).filter(
        Publications.doi.in_(doi_list)).all()
    matching_db_dois = [p.doi for p in matching_db_publications]

    new_pub_list = [p for p in pub_list if p.doi not in matching_db_dois]

    return new_pub_list, matching_db_publications


# merge 2 sets of publications
def merge_2_sets_of_publications(pub_set_1: List[Publications], pub_set_2: List[Publications]) -> \
        List[Publications]:
    pub_set_2_dois = [p.doi for p in pub_set_2]

    # check if a publication is listed in both the sets by comparing the doi
    common_publication_dois = [p.doi for p in pub_set_1 if p.doi in pub_set_2_dois]

    if len(common_publication_dois) == 0:
        final_pub_list = pub_set_1 + pub_set_2
    else:
        # remove the publications found in both sets from set 1
        updated_pub_set_1 = [p for p in pub_set_1 if p.doi not in common_publication_dois]
        final_pub_list = updated_pub_set_1 + pub_set_2

    return final_pub_list


# merge user provided publications with litvar publications. Merge all with the publications already stored in the
# database.
def merge_user_and_litvar_and_db_publications(user_pub_list: List[Publications],
                                              litvar_publications: List[Publications],
                                              db_publications: List[Publications]) -> \
        List[Publications]:
    # merge user provided publications with litvar publications
    user_and_litvar_pub_list = merge_2_sets_of_publications(user_pub_list, litvar_publications)

    # merge [user provided publication links & the litvar publications] and the publications already stored in the db
    final_pub_list = merge_2_sets_of_publications(user_and_litvar_pub_list, db_publications)

    return final_pub_list


#TODO: stop giving this function publications which the variant already has in first parameter
def store_variant_publications_in_db(publications: List[Publications], variant_publications: List[Publications]):
    # check which publications already exist in the publications table in the db
    vus_new_publications, vus_existing_publications = extract_publications_already_stored_in_db(publications)

    # store new publications
    for p in vus_new_publications:
        db.session.add(p)

        # set up relationship between the variant & its publications
        variant_publications.append(p)

    # set up relationships of variants with existing publications
    for p in vus_existing_publications:
        variant_publications.append(p)


def retrieve_and_store_variant_publications(vus_df: pd.DataFrame, variants_already_stored_in_db: bool):
    current_app.logger.info('Retrieving publications for variants')

    # iterate through the dataframe
    for index, row in vus_df.iterrows():
        publication_links: List[Publications] = []

        # check if user provided any literature links
        if 'Literature Links' in vus_df.keys():
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

        # attempt to get litvar publications using rsid if it exists, else use hgvs
        if row['RSID'] == 'NORSID':
            litvar_rsid = None
        else:
            litvar_rsid = row['RSID']

        # retrieve LitVar publications
        litvar_publications_res = get_publications(row['HGVS'], litvar_rsid, None)

        if litvar_publications_res.status != 200:
            current_app.logger.error(
                f'Retrieval of information for the user provided literature link failed 500')
            return InternalResponse(None, 500)
        else:
            litvar_publications = litvar_publications_res.data

        # retrieve variant
        variant = get_variant_from_db(row)

        # merge the user's publications together with litvar's publications together with the variant's db publications
        if variants_already_stored_in_db:
            final_pub_list = merge_user_and_litvar_and_db_publications(publication_links, litvar_publications,
                                                                       variant.publications)
        # merge the user's publications together with litvar's publications
        else:
            final_pub_list = merge_2_sets_of_publications(publication_links, litvar_publications)

        # extract the publications not yet included for the variant
        pub_not_in_variant = [p for p in final_pub_list if p not in variant.publications]
        store_variant_publications_in_db(pub_not_in_variant, variant.publications)

    return InternalResponse(None, 200)


def get_publications_by_variant_id_from_db(variant_id: str) -> (Variants, List[Dict]):
    variant: Variants = db.session.query(Variants).filter(Variants.id == variant_id).one_or_none()

    publication_list = []

    var_publications: List[Publications] = variant.publications

    for p in var_publications:
        encoded_publication = alchemy_encoder(p)

        # changing keys of dictionary
        encoded_publication['publicationId'] = encoded_publication.pop('id')
        date_published = encoded_publication.pop('date_published')
        encoded_publication['isSupplementaryMaterialMatch'] = encoded_publication.pop('match_in_sup_material')

        if date_published is not None:
            encoded_publication['date'] = date_published.strftime('%Y/%m/%d')

        if p.authors is not None:
            encoded_publication['authors'] = p.authors.replace('\"', '').strip('{}').split(',')

        publication_list.append(encoded_publication)

    return variant, publication_list


def check_for_new_litvar_publications():
    variants = db.session.query(Variants).all()

    for v in variants:
        # get rsid, if it exists
        db_snp_external_ref: ExternalReferences = db.session.query(ExternalReferences).filter(ExternalReferences.variant_id == v.id,
                                                                          ExternalReferences.db_type == 'db_snp').one_or_none()

        rsid = None
        if db_snp_external_ref is not None:
            rsid = db_snp_external_ref.db_snp.rsid

        # get one of its HGVS
        hgvs = v.variant_hgvs[0].hgvs

        # retrieve LitVar publications using rsid if it exists, else use hgvs
        litvar_publications_res = get_publications(hgvs, rsid, None)

        if litvar_publications_res.status != 200:
            current_app.logger.error(
                f'Retrieval of information for the user provided literature link failed 500')
        else:
            litvar_publications = litvar_publications_res.data

            # merge litvar publications with those found in db
            updated_pub_list = merge_2_sets_of_publications(litvar_publications, v.publications)

            # extract the publications not yet included for the variant
            pub_not_in_variant = [p for p in updated_pub_list if p not in v.publications]

            # update the variant's publications in the db
            store_variant_publications_in_db(pub_not_in_variant, v.publications)

    try:
        # Commit the session to persist changes to the database
        db.session.commit()
        return InternalResponse({'isSuccess': True}, 200)
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since insertion of new publications for variant with id {v.id} in DB failed '
            f'due to error: {e}')
        return InternalResponse({'isSuccess': False}, 500)
