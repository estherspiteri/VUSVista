from flask import Blueprint
import requests
from Bio import Entrez
import json
from typing import List

entrez_views = Blueprint('entrez_views', __name__)


@entrez_views.route('/pmids/<string:pmids>') #pmids seperated by commas
def retrieve_pubmed_publications_info(pmids: str):    
    Entrez.email = "esther.spiteri.18@um.edu.mt"
    publications_record = {}

    try:
        # retrieve unique identifiers of variants
        publications_handle = Entrez.efetch(db="pubmed", id={pmids})
        publications_record = Entrez.read(publications_handle)
        publications_handle.close()
    except IOError:
        print('Network error when calling Entrez.efetch()!')

    return publications_record