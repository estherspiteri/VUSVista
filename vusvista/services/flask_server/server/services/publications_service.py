import json
import math
import os
from datetime import datetime
import re
from typing import List, Dict
from urllib.parse import urlencode

import pandas as pd
import requests
from flask import current_app, Response, render_template
from flask_mail import Message
from sqlalchemy import desc
from sqlalchemy.exc import SQLAlchemyError

from server import db, mail
from server.helpers.data_helper import alchemy_encoder
from server.models import Publications, Variants, ExternalReferences, VariantsPublications, AutoPublicationEvalDates, \
    ScientificMembers
from server.responses.internal_response import InternalResponse
from server.services.litvar_service import get_publications


def get_publication_info(publication_link: str) -> InternalResponse:
    # check if there PMC-OA identifier
    pattern = r'PMC\d+'
    match = re.search(pattern, publication_link)
    if match:
        # extract the PMC identifier
        pmc_id = match.group()
        url_encoded_pmc_id = urlencode({'id': pmc_id})
        url = f"https://bg.api.oa.works/find?{url_encoded_pmc_id}"
    else:
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
            date = metadata['year'] + '-01-01'

        if date is not None and re.match(r"^\d{4}-\d{2}-\d{2}$", date):
            date = datetime.strptime(date, '%Y-%m-%d')
        else:
            date = None

        # eliminate any keywords included
        updated_doi = metadata.get('doi', None)
        if updated_doi is not None:
            updated_doi = updated_doi.replace("\"", "").split(',')[0]

        title = metadata.get('title', None)
        if title is not None:
            title = title[0].upper() + title[1:]

        publication = Publications(title=title, pmid=metadata.get('pmid', None),
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


def store_variant_publications_in_db(publications: List[Publications], variant_id: int,
                                     manually_added_pub_dois: List[str], date_added=datetime.now()):
    # check which publications already exist in the publications table in the db
    vus_new_publications, vus_existing_publications = extract_publications_already_stored_in_db(publications)

    # store new publications
    for p in vus_new_publications:
        db.session.add(p)
        db.session.flush()

    all_publications = vus_new_publications + vus_existing_publications

    # set up relationship between the variant & its publications
    for p in all_publications:
        vus_pub = VariantsPublications(variant_id=variant_id, publication_id=p.id, date_added=date_added,
                                       is_manually_added=p.doi is not None and p.doi in manually_added_pub_dois)
        db.session.add(vus_pub)


def retrieve_and_store_variant_publications(vus_df: pd.DataFrame):
    current_app.logger.info('Retrieving publications for variants')

    # iterate through the dataframe
    for index, row in vus_df.iterrows():
        publication_links: List[Publications] = []

        # check if user provided any literature links
        if 'Literature Links' in vus_df.keys():
            links = []

            if str(row['Literature Links']) != 'nan':
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

        # attempt to get litvar publications using rsid if it exists and if it does not have any errors, else use hgvs
        if row['RSID'] == 'NORSID' or len(row['RSID dbSNP errorMsgs']) > 0:
            litvar_rsid = None
        else:
            litvar_rsid = row['RSID']

        hgvs = None
        if not (isinstance(row["HGVS"], float) and math.isnan(row["HGVS"])) and row["HGVS"] is not None:
            hgvs = row['HGVS'].split(' ')[0]

        litvar_publications = []
        if hgvs is not None or litvar_rsid is not None:
            # retrieve LitVar publications
            litvar_publications_res = get_publications(hgvs, litvar_rsid, None)

            if litvar_publications_res.status != 200:
                current_app.logger.error(
                    f'Retrieval of information for the user provided literature link failed 500')
            else:
                litvar_publications = litvar_publications_res.data

        # retrieve variant
        variant: Variants = db.session.query(Variants).get(row['Variant Id'])

        # get variant's current publications
        variant_pub_ids = [vp.publication_id for vp in variant.variants_publications]
        variant_publications: List[Publications] = db.session.query(Publications).filter(
            Publications.id.in_(variant_pub_ids)).all()

        # merge the user's publications together with litvar's publications together with the variant's db publications
        final_pub_list = merge_user_and_litvar_and_db_publications(publication_links, litvar_publications,
                                                                       variant_publications)

        date = datetime.now()

        # add new evaluation date
        auto_pub_eval_dates = AutoPublicationEvalDates(variant_id=variant.id, eval_date=date)
        db.session.add(auto_pub_eval_dates)

        # get doi of manually added links
        manually_added_dois = [p.doi for p in publication_links]

        # extract the publications not yet included for the variant
        variant_publication_dois = [vp.publication.doi for vp in variant.variants_publications]
        pub_not_in_variant = [p for p in final_pub_list if p.doi not in variant_publication_dois]
        store_variant_publications_in_db(pub_not_in_variant, variant.id, manually_added_dois, date)

    return InternalResponse(None, 200)


def get_publications_by_variant_id_from_db(variant_id: str) -> (Variants, List[Dict]):
    variant: Variants = db.session.query(Variants).filter(Variants.id == variant_id).one_or_none()

    publication_list = []

    # get variant's current publications
    variant_pub_ids = [vp.publication_id for vp in variant.variants_publications]
    var_publications: List[Publications] = db.session.query(Publications).filter(
        Publications.id.in_(variant_pub_ids)).all()

    # get publications ids of those publications added manually to variant
    manually_added_pub_ids = [vp.publication_id for vp in variant.variants_publications if vp.is_manually_added]

    for p in var_publications:
        encoded_publication = alchemy_encoder(p)

        # changing keys of dictionary
        encoded_publication['publicationId'] = encoded_publication.pop('id')
        date_published = encoded_publication.pop('date_published')
        encoded_publication['isSupplementaryMaterialMatch'] = encoded_publication.pop('match_in_sup_material')
        encoded_publication['isAddedManually'] = p.id in manually_added_pub_ids

        if date_published is not None:
            encoded_publication['date'] = date_published.strftime('%Y/%m/%d')

        if p.authors is not None:
            encoded_publication['authors'] = p.authors.replace('\"', '').strip('{}').split(',')

        publication_list.append(encoded_publication)

    return variant, publication_list


def update_variant_publications(variant: Variants, hgvs: str | None, rsid: str | None) -> List:
    pub_not_in_variant = []
    if hgvs is not None or rsid is not None:
        # retrieve LitVar publications using rsid if it exists, else use hgvs
        litvar_publications_res = get_publications(hgvs, rsid, None)

        if litvar_publications_res.status != 200:
            current_app.logger.error(
                f'Retrieval of information for the user provided literature link failed 500')
        else:
            litvar_publications = litvar_publications_res.data

            if len(litvar_publications) > 0:
                # get variant's current publications
                variant_publications = [vp.publication for vp in variant.variants_publications]
                variant_publications_doi = [p.doi for p in variant_publications]

                # merge litvar publications with those found in db
                updated_pub_list = merge_2_sets_of_publications(litvar_publications, variant_publications)

                # extract the publications not yet included for the variant
                pub_not_in_variant = [p for p in updated_pub_list if p.doi not in variant_publications_doi]

                date = datetime.now()

                # update the variant's publications in the db
                store_variant_publications_in_db(pub_not_in_variant, variant.id, [], date)

            auto_pub_eval_date = AutoPublicationEvalDates(eval_date=date, variant_id=variant.id)
            db.session.add(auto_pub_eval_date)
    return pub_not_in_variant


def check_for_new_litvar_publications():
    variants = db.session.query(Variants).all()
    updates = []

    for v in variants:
        # get rsid, if it exists
        db_snp_external_ref: ExternalReferences = db.session.query(ExternalReferences).filter(
            ExternalReferences.variant_id == v.id,
            ExternalReferences.db_type == 'db_snp').one_or_none()

        rsid = None
        if db_snp_external_ref is not None and len(db_snp_external_ref.error_msg) == 0:
            rsid = db_snp_external_ref.db_snp.rsid

        hgvs = None
        # get one of its HGVS
        if len(v.variant_hgvs) > 0:
            hgvs = v.variant_hgvs[0].hgvs.split(' ')[0]

        pub_not_in_variant = update_variant_publications(v, hgvs, rsid)

        if len(pub_not_in_variant) > 0:
            pub_not_in_variant_info = [{"link": p.link, "doi": p.doi, "title": p.title} for p in pub_not_in_variant]
            update = ({'variant_id': v.id, 'chromosome': v.chromosome,
                       'chromosome_position': v.chromosome_position,
                       'gene': v.gene_name, 'alt_allele': v.alt,
                       'ref_allele': v.ref, "publications": pub_not_in_variant_info})

            updates.append(update)

    try:
        # Commit the session to persist changes to the database
        db.session.commit()

        if len(updates) > 0:
            # get all scientific members of staff
            recipients: List[ScientificMembers] = db.session.query(ScientificMembers).all()
            recipient_emails = [r.email for r in recipients]

            mail_message = Message(
                subject='Publication Updates',
                bcc=recipient_emails,
                html=render_template('publication_update_email_template.html', updates=updates,
                                     domain=os.environ.get('FRONT_URL')),
            )

            try:
                current_app.logger.info('Sending email with Publication updates.')
                mail.send(mail_message)
            except Exception as e:
                current_app.logger.error(f'Publication updates email was not sent successfully: {e}')

        return InternalResponse({'isSuccess': True}, 200)
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since insertion of new publications for variant with id {v.id} in DB failed '
            f'due to error: {e}')
        return InternalResponse({'isSuccess': False}, 500)


def get_variant_publication_updates(variant_id: str):
    variant: Variants = db.session.get(Variants, int(variant_id))

    variant_publication_entries: List[VariantsPublications] = variant.variants_publications

    auto_publication_updates = {}
    manual_publication_updates = {}
    for vp in variant_publication_entries:
        publication: Publications = db.session.query(Publications).get(vp.publication_id)
        date = datetime.strftime(vp.date_added, '%d/%m/%Y %H:%M')

        if vp.is_manually_added:
            update_arr = manual_publication_updates.get(date, [])
            manual_publication_updates[date] = update_arr + [
                {"title": publication.title, "doi": publication.doi, "link": publication.link, "isManuallyAdded": True}]
        else:
            update_arr = auto_publication_updates.get(date, [])
            auto_publication_updates[date] = update_arr + [
                {"title": publication.title, "doi": publication.doi, "link": publication.link,
                 "isManuallyAdded": False}]

    auto_pub_eval_dates: List[AutoPublicationEvalDates] = db.session.query(AutoPublicationEvalDates).filter(
        AutoPublicationEvalDates.variant_id == variant_id).all()

    # get both manual and automatic updates dates
    all_auto_pub_eval_dates = [d.eval_date for d in auto_pub_eval_dates]
    all_manual_pub_eval_dates = [vp.date_added for vp in variant_publication_entries if vp.is_manually_added == True]
    all_dates = all_auto_pub_eval_dates + all_manual_pub_eval_dates

    # reversed to get dates in desc order
    all_dates.sort()
    all_dates.reverse()

    all_dates_str = [datetime.strftime(d, '%d/%m/%Y %H:%M') for d in all_dates]
    all_dates_str_unique = []
    for d in all_dates_str:
        if d not in all_dates_str_unique:
            all_dates_str_unique.append(d)

    all_auto_pub_eval_dates_str = [datetime.strftime(d.eval_date, '%d/%m/%Y %H:%M') for d in auto_pub_eval_dates]

    variant_publications_update_list = []
    for date in all_dates_str_unique:
        if date in all_auto_pub_eval_dates_str and date not in auto_publication_updates.keys():
            pub_updates = None

        else:
            pub_updates = []

            if date in auto_publication_updates.keys():
                pub_updates += auto_publication_updates[date]

            if date in manual_publication_updates.keys():
                pub_updates += manual_publication_updates[date]

        variant_publications_update_list.append({'lastEval': date,
                                                 "publicationUpdates": pub_updates})

    dates_with_updates = list(set([u['lastEval'].split(" ")[0] for u in variant_publications_update_list if
                                   u['publicationUpdates'] is not None]))

    return Response(json.dumps({'isSuccess': True, 'variantPublicationUpdates': variant_publications_update_list,
                                'datesWithUpdates': dates_with_updates}), 200, mimetype='application/json')


def add_publications_to_variant(variant_id: str, publication_links_list: List[str]):
    variant: Variants = db.session.query(Variants).get(variant_id)
    variant_publication_ids = [vp.publication_id for vp in variant.variants_publications]
    variant_publications: List[Publications] = db.session.query(Publications).filter(
        Publications.id.in_(variant_publication_ids)).all()

    added_publications: List[Publications] = []
    for link in publication_links_list:
        link_publication_res = get_publication_info(link)

        if link_publication_res.status != 200:
            current_app.logger.error(
                f'Retrieval of information for the user provided literature link failed 500')
        else:
            new_pub = link_publication_res.data
            added_publications.append(new_pub)

    final_pub_list = merge_2_sets_of_publications(added_publications, variant_publications)

    # extract the publications not yet included for the variant
    pub_not_in_variant = [p for p in final_pub_list if p.id not in variant_publication_ids]
    pub_not_in_variant_dois = [p.doi for p in pub_not_in_variant]
    store_variant_publications_in_db(pub_not_in_variant, variant.id, pub_not_in_variant_dois)

    try:
        # Commit the session to persist changes to the database
        db.session.commit()

        variant, publications = get_publications_by_variant_id_from_db(variant_id)
        return Response(json.dumps({'isSuccess': True, 'publications': publications}), 200, mimetype='application/json')
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since insertion of new publication entries in DB failed due to error: {e}')
        return Response(json.dumps({'isSuccess': False, "publications": None}), 200, mimetype='application/json')


def remove_publications_to_variant(variant_id: str, publication_id: int):
    delete_publication = False

    variant_publication: VariantsPublications = db.session.query(VariantsPublications).filter(VariantsPublications.variant_id == variant_id, VariantsPublications.publication_id == publication_id).one()

    # delete publication if it only has the variant that it will be deleted from as a variant
    if len(variant_publication.publication.variants_publications) == 1:
        delete_publication = True

    db.session.delete(variant_publication)

    if delete_publication:
        db.session.query(Publications).filter(Publications.id == publication_id).delete()

    try:
        # Commit the session to persist changes to the database
        db.session.commit()

        variant, publications = get_publications_by_variant_id_from_db(variant_id)
        return Response(json.dumps({'isSuccess': True, 'publications': publications}), 200, mimetype='application/json')
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since removal of publication {publication_id} from variant {variant_id} in DB failed due to error: {e}')
        return Response(json.dumps({'isSuccess': False, "publications": None}), 200, mimetype='application/json')


def get_last_pub_auto_update_date() -> str | None:
    auto_pub_eval_date: AutoPublicationEvalDates = db.session.query(AutoPublicationEvalDates).order_by(
        desc(AutoPublicationEvalDates.eval_date)).first()

    if auto_pub_eval_date is not None:
        return datetime.strftime(auto_pub_eval_date.eval_date, '%d/%m/%Y %H:%M')
    else:
        return None
