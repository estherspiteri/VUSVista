from typing import Dict
import pandas as pd

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