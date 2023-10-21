from flask import Blueprint, url_for, request
import requests
import pandas as pd
from urllib.parse import urlencode
from typing import Dict
from helpers import extract_abstracts_by_pmids, add_abstracts_to_df


views = Blueprint('views', __name__)

@views.route('/publications/<string:rsid>')
def get_publications_of_variant(rsid: str):
    domain = request.host

    # get LitVar id
    litvar_id_url = url_for('litvar_views.get_litvar_id', rsid=rsid)
    litvar_id_res = requests.get(f'http://{domain}{litvar_id_url}')

    if litvar_id_res.status_code != 200:
        print(
            f'LitVar Search Variant query failed! - {litvar_id_res.reason}')
        exit(1)
    else:
        litvar_id = litvar_id_res.text

        if litvar_id is not None:
            # search for LitVar publications for given variant
            litvar_publications_url = url_for('litvar_views.get_litvar_publications', litvar_id=litvar_id)
            litvar_publications_res = requests.get(f'http://{domain}{litvar_publications_url}')

            if litvar_publications_res.status_code != 200:
                print(
                f'LitVar Variant Publications query failed! - {litvar_publications_res.reason}')
                exit(1)
            else:
                litvar_publications = litvar_publications_res.json()
                publications_df = pd.DataFrame.from_records(litvar_publications['results'])

                # extract pmids
                litvar_pmids = [str(id) for id in publications_df['pmid']]

                # retrieve more information about publications
                pubmed_publications_url = url_for('entrez_views.retrieve_pubmed_publications_info', pmids=','.join(litvar_pmids))
                pubmed_publications_res = requests.get(f'http://{domain}{pubmed_publications_url}')
                
                if pubmed_publications_res.status_code != 200:
                    print(
                    f'Entrez Publications query failed! - {pubmed_publications_res.reason}')
                    exit(1)
                else:
                    pubmed_publications_info = pubmed_publications_res.json()

                    # store publication abstracts in a dict where the key is the publication's PMID
                    pubmed_publications_abstracts_dict = extract_abstracts_by_pmids(pubmed_publications_info)

                    # add abstract column to LitVar publications df
                    publications_df = add_abstracts_to_df(publications_df, pubmed_publications_abstracts_dict)

                    return publications_df.to_dict()
        else:
            return None
