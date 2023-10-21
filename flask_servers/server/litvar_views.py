from flask import Blueprint
import requests
from urllib.parse import urlencode
import pandas as pd


litvar_views = Blueprint('litvar_views', __name__)


@litvar_views.route('/id/<string:rsid>')
def get_litvar_id(rsid: str) -> str | None:
    litvar_search_variant_res = requests.get(
        f'https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/autocomplete/?query={rsid}')

    if litvar_search_variant_res.status_code != 200:
        print(
            f'LitVar Search Variant query failed! - {litvar_search_variant_res.reason}')
        exit(1)
    else:
        litvar_search_variant_res_json = litvar_search_variant_res.json()

        if len(litvar_search_variant_res_json) > 0:
            # TODO: check gene match

            # assuming first result is the most relevant
            return litvar_search_variant_res_json[0]['_id']
        else:
            print(
                f'LitVar Search Variant query - no LitVar id found for RSID {rsid}!')
            return None

@litvar_views.route('/publications/<path:litvar_id>')
def get_litvar_publications(litvar_id: str):# TODO: litvar_id needs to be url encoded before hand for the hash symbols to be recognized
    url_encoded_id = urlencode({'query': litvar_id}).split('=')[1]

    url = f"https://www.ncbi.nlm.nih.gov/research/litvar2-api/search/?variant={url_encoded_id}&sort=score%20desc"
    print(litvar_id)
    litvar_publications_res = requests.get(url)

    if litvar_publications_res.status_code != 200:
        print(f'LitVar Variant Publications query failed! - {litvar_publications_res.reason}')
        exit(1)
    else:
        return litvar_publications_res.json()
