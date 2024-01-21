from flask import current_app, Response
import requests
import pandas as pd
from urllib.parse import urlencode

from server.helpers.data_helper import convert_df_to_list
from server.responses.internal_response import InternalResponse
from server.services.entrez_service import retrieve_pubmed_publications_info
from typing import Dict
import json


def get_litvar_id(rsid: str) -> InternalResponse:
    litvar_search_variant_res = requests.get(
        f'https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/autocomplete/?query={rsid}')

    if litvar_search_variant_res.status_code != 200:
        current_app.logger.error(
            f'Response failure {litvar_search_variant_res.status_code}: {litvar_search_variant_res.reason}')
        return InternalResponse(None, litvar_search_variant_res.status_code, litvar_search_variant_res.reason)
    else:
        litvar_search_variant_res_json = litvar_search_variant_res.json()

        if len(litvar_search_variant_res_json) > 0:
            # TODO: check gene match

            # assuming first result is the most relevant
            litvar_id = litvar_search_variant_res_json[0]['_id']

            current_app.logger.info(f"Litvar ID for variant rsid {rsid} is {litvar_id}")
            return InternalResponse(litvar_id, 200)
        else:
            current_app.logger.info(f'LitVar Search Variant query - no LitVar id found for RSID {rsid}!')
            return InternalResponse('', litvar_search_variant_res.status_code, litvar_search_variant_res.reason)


def get_litvar_publications(litvar_id: str) -> InternalResponse:
    url_encoded_id = urlencode({'query': litvar_id}).split('=')[1]
    url = f"https://www.ncbi.nlm.nih.gov/research/litvar2-api/search/?variant={url_encoded_id}&sort=score%20desc"

    litvar_publications_res = requests.get(url)

    if litvar_publications_res.status_code != 200:
        current_app.logger.error(
            f'Response failure {litvar_publications_res.status_code}: {litvar_publications_res.reason}')
        return InternalResponse(None, litvar_publications_res.status_code, litvar_publications_res.reason)
    else:
        return InternalResponse(litvar_publications_res.json(), 200)


def extract_abstracts_by_pmids(pubmed_publications_info) -> Dict:
    abstract_dict = {}

    for pubmed_article_info in pubmed_publications_info['PubmedArticle']:
        # extract pmid
        pmid = int(pubmed_article_info['MedlineCitation']['PMID'])

        if 'Abstract' in pubmed_article_info['MedlineCitation']['Article']:
            # extract abstract and concatenate parts of abstract using '\n'
            abstract = '\n'.join([str(abstract_line) for abstract_line in pubmed_article_info['MedlineCitation']['Article']['Abstract']['AbstractText']])
        else:
            abstract = ''

        # set key-value pair
        abstract_dict[pmid] = abstract

    return abstract_dict


def add_abstracts_to_df(publications_df: pd.DataFrame, abstract_dict: Dict):
    for pmid in abstract_dict.keys():
        # extract rows that match the given pmid using boolean indexing
        matching_rows = publications_df['pmid'] == pmid

        # add the abstract to the rows
        publications_df.loc[matching_rows, 'abstract'] = abstract_dict[pmid]

    return publications_df


def get_publications(rsid: str) -> Response:
    # get LitVar id
    current_app.logger.info(f'Retrieving LitVar ID for RSID {rsid}')
    litvar_id_res: InternalResponse = get_litvar_id(rsid)

    if litvar_id_res.status != 200:
        current_app.logger.error(
            f'LitVar Search Variant query failed 500')
        return Response(json.dumps({'isSuccess': False}), 500)
    else:
        litvar_id = litvar_id_res.data

        if len(litvar_id) > 0:
            current_app.logger.info(f'Retrieved LitVar ID {litvar_id}')

            # search for LitVar publications for given variant
            current_app.logger.info(f'Retrieving LitVar publications for {litvar_id}')
            litvar_publications_res: InternalResponse = get_litvar_publications(litvar_id)

            if litvar_publications_res.status != 200:
                current_app.logger.error(
                    f'LitVar Variant Publications query failed 500')
                return Response(json.dumps({'isSuccess': False}), 500)
            else:
                litvar_publications = litvar_publications_res.data
                publications_df = pd.DataFrame.from_records(litvar_publications['results'])

                # extract pmids
                litvar_pmids = [str(id) for id in publications_df['pmid']]

                # retrieve more information about publications
                current_app.logger.info(f'Retrieving PubMed information for LitVar publications')
                pubmed_publications_res: InternalResponse = retrieve_pubmed_publications_info(','.join(litvar_pmids))

                if pubmed_publications_res.status != 200:
                    current_app.logger.error(
                        f'Entrez Publications query failed 500')
                    return Response(json.dumps({'isSuccess': False}), 500)
                else:
                    current_app.logger.info(f'Appending PubMed abstracts to LitVar publications')

                    pubmed_publications_info = pubmed_publications_res.data

                    # store publication abstracts in a dict where the key is the publication's PMID
                    pubmed_publications_abstracts_dict = extract_abstracts_by_pmids(pubmed_publications_info)

                    # add abstract column to LitVar publications df
                    publications_df = add_abstracts_to_df(publications_df, pubmed_publications_abstracts_dict)

                    # replace NaNs with empty strings
                    publications_df = publications_df.fillna('')

                    publications_list = convert_df_to_list(publications_df)

                    current_app.logger.info(f'Sending user {len(publications_list)} publications')
                    return Response(json.dumps({'isSuccess': True, 'publicationSearch': {'publications': publications_list, "isLitvarIdFound": True}}), 200)
        else:
            current_app.logger.info(f'Sending user 0 publications')
            return Response(json.dumps({'isSuccess': True, 'publicationSearch': {'publications': [], "isLitvarIdFound": False}}), 200)
