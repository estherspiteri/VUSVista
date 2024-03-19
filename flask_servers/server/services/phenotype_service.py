from typing import Dict

import requests
from flask import current_app
from sqlalchemy.exc import SQLAlchemyError

from server import db
from server.models import Phenotypes, Samples
from server.responses.internal_response import InternalResponse
import urllib
from urllib.parse import urlencode


def get_hpo_terms(phenotype: str) -> InternalResponse:
    url_encoded_phenotype = urllib.parse.quote(phenotype)
    url = f"https://hpo.jax.org/api/hpo/search?q={url_encoded_phenotype}&max=100&category=terms"

    hpo_res = requests.get(url)

    if hpo_res.status_code != 200:
        current_app.logger.error(
            f'Response failure {hpo_res.status_code}: {hpo_res.reason}')
        return InternalResponse({'isSuccess': False, 'hpoTerms': None}, 500)
    else:
        terms = [{'ontologyId': x['ontologyId'], 'name': x['name']} for x in hpo_res.json()['terms']]
        return InternalResponse({'isSuccess': True, 'hpoTerms': terms}, 200)


def get_hpo_term_from_phenotype_name(name: str) -> InternalResponse:
    terms_res = get_hpo_terms(name)

    if terms_res.status != 200:
        return terms_res
    else:
        term_arr = terms_res.data['hpoTerms']

        if len(term_arr) == 1:
            current_app.logger.info(f'HPO term found for phenotype name: {name}')
            return InternalResponse(term_arr[0], 200)
        elif len(term_arr) > 1:
            terms_filtered_arr = [t for t in term_arr if t['name'].lower() == name.lower()]

            if len(terms_filtered_arr) == 1:
                return InternalResponse(terms_filtered_arr[0], 200)
            elif len(terms_filtered_arr) > 1:
                current_app.logger.error(f'Too many hpo terms found for  phenotype name: {name}')

        current_app.logger.error(f'No matching hpo terms for phenotype name: {name}')
        return InternalResponse({}, 500)


def append_phenotype_to_sample(sample: Samples, phenotype: Dict):
    phenotype_selected: Phenotypes = db.session.query(Phenotypes).filter(
        Phenotypes.ontology_term_id == phenotype['ontologyId'],
    ).one_or_none()

    if phenotype_selected is None:
        new_phenotype = Phenotypes(ontology_term_id=phenotype['ontologyId'], term_name=phenotype['name'])
        db.session.add(new_phenotype)

        sample.ontology_term.append(new_phenotype)
    else:
        sample.ontology_term.append(phenotype_selected)


def add_phenotype_to_existing_sample(sample_id: str, phenotype: Dict):
    sample: Samples = db.session.query(Samples).filter(Samples.id == sample_id).first()

    append_phenotype_to_sample(sample, phenotype)

    try:
        # Commit the session to persist changes to the database
        db.session.commit()
        return InternalResponse({'isSuccess': True}, 200)
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since insertion of sample phenotype in DB failed due to error: {e}')
        return InternalResponse({'isSuccess': False}, 500)


def remove_phenotype_to_sample(sample_id: str, phenotype: Dict):
    sample: Samples = db.session.query(Samples).filter(Samples.id == sample_id).first()

    phenotype_removed: Phenotypes = db.session.query(Phenotypes).filter(
        Phenotypes.ontology_term_id == phenotype['ontologyId'],
    ).first()

    sample.ontology_term.remove(phenotype_removed)

    try:
        # Commit the session to persist changes to the database
        db.session.commit()
        return InternalResponse({'isSuccess': True}, 200)
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since deletion of sample phenotype in DB failed due to error: {e}')
        return InternalResponse({'isSuccess': False}, 500)
