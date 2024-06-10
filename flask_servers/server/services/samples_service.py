from typing import List, Dict

from flask import current_app
from sqlalchemy.exc import SQLAlchemyError

from server import db
from server.helpers.data_helper import get_variant_summary
from server.models import Samples, VariantsSamples, t_samples_phenotypes, Phenotypes, \
    Variants, VariantsSamplesUploads, FileUploads
from server.responses.internal_response import InternalResponse


def get_sample_info_from_db(sample: Samples) -> Dict:
    # Query with filter condition on sample_id
    phenotype_ontology_term_ids_res: List[str] = (db.session.query(t_samples_phenotypes.c.ontology_term_id)
                                                  .filter(t_samples_phenotypes.c.sample_id == sample.id).all())

    phenotype_ontology_term_ids = [x[0] for x in phenotype_ontology_term_ids_res]

    sample_phenotypes: List[Phenotypes] = (db.session.query(Phenotypes)
                                           .filter(Phenotypes.ontology_term_id.in_(phenotype_ontology_term_ids)).all())

    phenotypes = []

    for p in sample_phenotypes:
        sample_phenotype = {'ontologyId': p.ontology_term_id, 'name': p.term_name}
        phenotypes.append(sample_phenotype)

    variants_samples: List[VariantsSamples] = db.session.query(VariantsSamples).filter(
        VariantsSamples.sample_id == sample.id).all()

    variants = []

    for v_s in variants_samples:
        variant_details: Variants = db.session.query(Variants).filter(Variants.id == v_s.variant_id).first()

        variant_summary = get_variant_summary(variant_details)

        variant_sample = {'variantId': v_s.variant_id, 'variant': variant_summary, 'genotype': v_s.genotype.value,
                          'hgvs': v_s.variant_hgvs.hgvs}
        variants.append(variant_sample)

    return {'sampleId': sample.id, 'phenotype': phenotypes, 'genomeVersion': sample.genome_version,
            'variants': variants}


def retrieve_all_samples_from_db():
    samples_arr = []

    # retrieve all samples
    all_samples = db.session.query(Samples).all()

    # for each sample retrieve the number of variants it has
    for sample in all_samples:
        variants_count: List[VariantsSamples] = db.session.query(VariantsSamples).filter(
            VariantsSamples.sample_id == sample.id).count()

        samples_arr.append({'sampleId': sample.id, 'numOfVariants': variants_count})

    return samples_arr


def retrieve_sample_from_db(sample_id: str):
    # retrieve sample
    sample = db.session.query(Samples).filter(Samples.id == sample_id).one_or_none()

    if sample is None:
        return InternalResponse({'isSuccess': True, 'sample_dict': None}, 404)

    sample_dict = get_sample_info_from_db(sample)

    return InternalResponse({'isSuccess': True, 'sample_dict': sample_dict}, 200)


def delete_sample_entry(sample_id: str) -> InternalResponse:
    sample: Samples = db.session.query(Samples).get(sample_id)

    if not sample:
        current_app.logger.error(f'Sample with ID: {sample_id} not found')
        return InternalResponse({'isSuccess': False}, 404)

    current_app.logger.info(
        f"Deleting sample with ID: {sample.id}")
    db.session.delete(sample)

    # deleting any file_uploads without variants_samples_uploads
    file_uploads_query = db.session.query(FileUploads).filter(~FileUploads.variants_samples_uploads.any())
    file_uploads_query.delete()

    # deleting any phenotypes without samples
    phenotypes_query = db.session.query(Phenotypes).filter(~Phenotypes.sample.any())
    phenotypes_query.delete()

    try:
        # Commit the session to persist changes to the database
        db.session.commit()
        return InternalResponse({'isSuccess': True}, 200)
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since deletion of sample {sample_id} in DB failed due to error: {e}')
        return InternalResponse({'isSuccess': False}, 500)
