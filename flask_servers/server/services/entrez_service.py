from Bio import Entrez
from flask import current_app

from server.responses.entrez_response import EntrezResponse


def retrieve_pubmed_publications_info(pmids: str):
    Entrez.email = "esther.spiteri.18@um.edu.mt"

    try:
        # retrieve unique identifiers of variants
        publications_handle = Entrez.efetch(db="pubmed", id={pmids})
        publications_record = Entrez.read(publications_handle)
        publications_handle.close()
    except IOError as e:
        current_app.logger.error(f"Network error when calling Entrez.efetch(): {str(e)}")
        return EntrezResponse(None, e.errno, str(e))

    return EntrezResponse(publications_record, 200, '')
