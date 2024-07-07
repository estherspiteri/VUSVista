from datetime import datetime
from urllib import parse

from flask import current_app
import requests

from server.models import Publications
from server.responses.internal_response import InternalResponse
from server.services.entrez_service import retrieve_pubmed_publications_info
from typing import Dict, List


def get_litvar_info(search_string: str) -> InternalResponse:
    litvar_search_variant_res = requests.get(
        f'https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/autocomplete/?query={search_string}')

    if litvar_search_variant_res.status_code != 200:
        current_app.logger.error(
            f'Response failure {litvar_search_variant_res.status_code}: {litvar_search_variant_res.reason}')
        return InternalResponse(None, 500, litvar_search_variant_res.reason)
    else:
        litvar_search_variant_res_json = litvar_search_variant_res.json()

        if len(litvar_search_variant_res_json) > 0:
            # TODO: check gene match

            # assuming first result is the most relevant
            litvar_id = litvar_search_variant_res_json[0]['_id']

            current_app.logger.info(f"Litvar ID for variant {search_string} is {litvar_id}")
            return InternalResponse(litvar_id, 200)
        else:
            current_app.logger.info(f'LitVar Search Variant query - no LitVar id found for {search_string}!')
            return InternalResponse('', 500, litvar_search_variant_res.reason)


# search by rsid if it exists else by hgvs (if it exists)
def get_litvar_id(hgvs: str | None, rsid: str | None) -> InternalResponse:
    formatted_hgvs = None
    if hgvs is not None:
        formatted_hgvs = hgvs.split(' ')[0]

    if rsid is not None:
        rsid_litvar_info_res = get_litvar_info(rsid)

        if rsid_litvar_info_res.status != 200:
            current_app.logger.error(
                f'No Litvar info found for RSID {rsid}')

            if hgvs is not None:
                hgvs_litvar_info_res = get_litvar_info(formatted_hgvs)

                if hgvs_litvar_info_res.status != 200:
                    current_app.logger.error(
                        f'No Litvar info found for HGVS {formatted_hgvs}')
                else:
                    current_app.logger.info(
                        f'Litvar info found for HGVS {formatted_hgvs}')
                return hgvs_litvar_info_res
        else:
            current_app.logger.info(
                f'Litvar info found for RSID {rsid}')
            return rsid_litvar_info_res
    else:
        hgvs_litvar_info_res = get_litvar_info(formatted_hgvs)

        if hgvs_litvar_info_res.status != 200:
            current_app.logger.error(
                f'No Litvar info found for HGVS {formatted_hgvs}')
        else:
            current_app.logger.info(
                f'Litvar info found for HGVS {formatted_hgvs}')
        return hgvs_litvar_info_res


def get_litvar_publications(litvar_id: str, optional_text: str | None) -> InternalResponse:
    url_encoded_id = parse.quote(litvar_id)

    if optional_text is None:
        url = f"https://www.ncbi.nlm.nih.gov/research/litvar2-api/search/?variant={url_encoded_id}&sort=score%20desc"
    else:
        url_encoded_optional_text = parse.quote(optional_text)
        url = f"https://www.ncbi.nlm.nih.gov/research/litvar2-api/search/?text={url_encoded_optional_text}&variant={url_encoded_id}&sort=score%20desc"

    litvar_publications_res = requests.get(url)

    if litvar_publications_res.status_code != 200:
        current_app.logger.error(
            f'Response failure {litvar_publications_res.status_code}: {litvar_publications_res.reason}')
        return InternalResponse(None, litvar_publications_res.status_code, litvar_publications_res.reason)
    else:
        litvar_publications_data = litvar_publications_res.json()

        litvar_publications: List[Publications] = []

        for publication in litvar_publications_data['results']: # TODO: include doi & fix date when received in front (convert in front)
            date = datetime.strptime(publication['date'], '%Y-%m-%dT%H:%M:%SZ').strftime('%Y/%m/%d')
            litvar_publications.append(Publications(title=publication['title'], pmid=publication['pmid'],
                                                    authors=publication['authors'], journal=publication['journal'],
                                                    date_published=date, match_in_sup_material=publication['is_sup_mat_match'],
                                                    link=f"https://pubmed.ncbi.nlm.nih.gov/{publication['pmid']}"))

        return InternalResponse(litvar_publications, 200)


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


def extract_doi_by_pmids(pubmed_publications_info) -> Dict:
    doi_dict = {}

    for pubmed_article_info in pubmed_publications_info['PubmedArticle']:
        # extract pmid
        pmid = int(pubmed_article_info['MedlineCitation']['PMID'])

        # extract doi
        id_list = pubmed_article_info['PubmedData']['ArticleIdList']
        doi = [x for x in id_list if x.startswith('10.')][0]

        # set key-value pair
        doi_dict[pmid] = doi

    return doi_dict


def add_abstracts_to_publications(publications: List[Publications], abstract_dict: Dict) -> List[Publications]:
    for publication in publications:
        # retrieve abstract based on publication's pmid
        publication.abstract = abstract_dict.get(publication.pmid, None)

    return publications


def add_doi_to_publications(publications: List[Publications], doi_dict: Dict) -> List[Publications]:
    for publication in publications:
        # retrieve doi based on publication's pmid
        publication.doi = doi_dict.get(publication.pmid, None)

    return publications


def get_more_info_on_litvar_publications(litvar_publications: list[Publications]):
    # extract pmids
    litvar_pmids = [str(p.pmid) for p in litvar_publications]

    # retrieve more information about publications
    current_app.logger.info(f'Retrieving PubMed information for LitVar publications')
    pubmed_publications_res: InternalResponse = retrieve_pubmed_publications_info(','.join(litvar_pmids))

    if pubmed_publications_res.status != 200:
        current_app.logger.error(
            f'Entrez Publications query failed 500')
        return InternalResponse({'isSuccess': False}, 500)
    else:
        current_app.logger.info(f'Appending PubMed abstracts to LitVar publications')

        pubmed_publications_info = pubmed_publications_res.data

        # store publication abstracts in a dict where the key is the publication's PMID
        pubmed_publications_abstracts_dict = extract_abstracts_by_pmids(pubmed_publications_info)

        # add abstract column to LitVar publications df
        litvar_publications = add_abstracts_to_publications(litvar_publications,
                                                            pubmed_publications_abstracts_dict)

        # store publication doi in a dict where the key is the publication's PMID
        pubmed_publications_doi_dict = extract_doi_by_pmids(pubmed_publications_info)

        # add abstract column to LitVar publications df
        litvar_publications = add_doi_to_publications(litvar_publications, pubmed_publications_doi_dict)

        current_app.logger.info(f'Found {len(litvar_publications)} Litvar publications')
        return InternalResponse(litvar_publications, 200)


def get_publications(hgvs: str | None, rsid: str | None, optional_text: str | None) -> InternalResponse:
    # get LitVar id
    current_app.logger.info(f'Retrieving LitVar ID for HGVS {hgvs} and RSID {rsid}')

    if hgvs is not None or rsid is not None:
        litvar_id_res: InternalResponse = get_litvar_id(hgvs, rsid)

        if litvar_id_res is None or litvar_id_res.status != 200:
            current_app.logger.error(
                f'LitVar Search Variant query failed 500')
            return InternalResponse({'isSuccess': False}, 500)
        else:
            litvar_id = litvar_id_res.data

            if len(litvar_id) > 0:
                current_app.logger.info(f'Retrieved LitVar ID {litvar_id}')

                # search for LitVar publications for given variant
                current_app.logger.info(f'Retrieving LitVar publications for {litvar_id}')
                litvar_publications_res: InternalResponse = get_litvar_publications(litvar_id, optional_text)

                if litvar_publications_res.status != 200:
                    current_app.logger.error(
                        f'LitVar Variant Publications query failed 500')
                    return InternalResponse({'isSuccess': False}, 500)
                else:
                    publications = litvar_publications_res.data

                    if len(publications) > 0:
                        # get additional information about the publications
                        handle_litvar_publications_res = get_more_info_on_litvar_publications(publications)

                        if handle_litvar_publications_res.status != 200:
                            current_app.logger.error(
                                f'LitVar Variant Publications query for more information failed 500')
                            return InternalResponse({'isSuccess': False}, 500)
                        else:
                            return InternalResponse(handle_litvar_publications_res.data, 200)

                current_app.logger.info(f'Found 0 LitVar publications')
    return InternalResponse([], 200)
