from Bio import Entrez
from flask import current_app

from server.responses.internal_response import InternalResponse


def retrieve_pubmed_publications_info(pmids: str) -> InternalResponse:
    Entrez.email = "esther.spiteri.18@um.edu.mt"

    try:
        # retrieve unique identifiers of variants
        publications_handle = Entrez.efetch(db="pubmed", id={pmids})
        publications_record = Entrez.read(publications_handle)
        publications_handle.close()
    except IOError as e:
        current_app.logger.error(f"Network error when calling Entrez.efetch(): {str(e)}")
        return InternalResponse(None, e.errno, str(e))

    return InternalResponse(publications_record, 200, '')
