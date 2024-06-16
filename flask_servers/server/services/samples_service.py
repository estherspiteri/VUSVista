from typing import List, Dict, Tuple

from flask import current_app
from sqlalchemy.exc import SQLAlchemyError

from server import db
from server.helpers.data_helper import get_variant_summary
from server.models import Samples, VariantsSamples, t_samples_phenotypes, Phenotypes, \
    Variants, FileUploads, VariantHgvs
from server.responses.internal_response import InternalResponse
from server.services.consequence_service import get_consequences_for_new_vus
from server.services.phenotype_service import append_phenotype_to_sample
from server.services.variants_samples_service import store_upload_details_for_variant_sample


def get_sample_variants(variants_samples: List[VariantsSamples]) -> Tuple[List, List]:
    variants = []

    for v_s in variants_samples:
        variant_details: Variants = db.session.query(Variants).filter(Variants.id == v_s.variant_id).first()

        variant_summary = get_variant_summary(variant_details)

        variant_sample = {'variantId': v_s.variant_id, 'variant': variant_summary, 'genotype': v_s.genotype.value,
                          'hgvs': v_s.variant_hgvs.hgvs, 'isHgvsUpdated': v_s.variant_hgvs.is_updated}
        variants.append(variant_sample)

    # get the variant that the sample does not have
    sample_variants_ids = [v.variant_id for v in variants_samples]

    not_sample_variants: List[Variants] = db.session.query(Variants).filter(
        Variants.id.not_in(sample_variants_ids)).all()

    not_sample_variants_list = []

    for v in not_sample_variants:
        variant_summary = get_variant_summary(v, True)

        variant_sample = {'variantId': v.id, 'variant': variant_summary}
        not_sample_variants_list.append(variant_sample)

    return variants, not_sample_variants_list


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

    variants, not_sample_variants_list = get_sample_variants(variants_samples)

    return {'sampleId': sample.id, 'phenotype': phenotypes, 'genomeVersion': sample.genome_version,
            'variants': variants, 'notSampleVariants': not_sample_variants_list}


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


def update_variant_sample_hgvs(sample_id: str, variant_id: str, hgvs: str):
    variant_sample: VariantsSamples = db.session.query(VariantsSamples).filter(VariantsSamples.sample_id == sample_id, VariantsSamples.variant_id == variant_id).first()
    variant_sample.variant_hgvs.hgvs = hgvs
    variant_sample.variant_hgvs.is_updated = True

    try:
        # Commit the session to persist changes to the database
        db.session.commit()
        return InternalResponse({'isSuccess': True}, 200)
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since update of hgvs {hgvs} in DB failed due to error: {e}')
        return InternalResponse({'isSuccess': False}, 500)


def add_variants_to_sample(sample_id: str, variants_to_add: List) -> InternalResponse:
    for v in variants_to_add:
        hgvs: VariantHgvs = db.session.query(VariantHgvs).filter(VariantHgvs.variant_id == v['variantId'], VariantHgvs.hgvs == v['hgvs']).one_or_none()

        if hgvs is None:
            hgvs = VariantHgvs(variant_id=v['variantId'], hgvs=v['hgvs'], is_updated=False)
            db.session.add(hgvs)
            db.session.flush()

        hgvs_consequence_dict = {}

        # get variant consequences though HGVS
        get_consequences_for_new_vus_res = get_consequences_for_new_vus([hgvs.hgvs])

        if get_consequences_for_new_vus_res.status != 200:
            current_app.logger.error(
                f'Retrieval of variant consequence failed 500')
        else:
            hgvs_consequence_dict = get_consequences_for_new_vus_res.data['consequences_dict']

        variant_sample = VariantsSamples(variant_id=v['variantId'], sample_id=sample_id, genotype=v['genotype'].upper(), variant_hgvs_id=hgvs.id, consequence=hgvs_consequence_dict.get(hgvs.hgvs, ""))

        db.session.add(variant_sample)

        store_upload_details_for_variant_sample(None, False, sample_id, v['variantId'])

    try:
        # Commit the session to persist changes to the database
        db.session.commit()

        # get updated sample's variants
        sample: Samples = db.session.query(Samples).get(sample_id)
        variants_samples: List[VariantsSamples] = sample.variants_samples
        updated_variants, updated_not_sample_variants = get_sample_variants(variants_samples)

        return InternalResponse({'isSuccess': True, 'updatedVariants': updated_variants, "updatedNotSampleVariants": updated_not_sample_variants}, 200)
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since addition of new variants to sample {sample_id} in DB failed due to error: {e}')
        return InternalResponse({'isSuccess': False}, 500)


def remove_variants_to_sample(sample_id: str, variant_ids_to_remove: List[str]) -> InternalResponse:
    db.session.query(VariantsSamples).filter(VariantsSamples.sample_id == sample_id,
                                                              VariantsSamples.variant_id.in_(
                                                                  variant_ids_to_remove)).delete()

    # deleting any file_uploads without variants_samples_uploads
    db.session.query(FileUploads).filter(~FileUploads.variants_samples_uploads.any()).delete()

    try:
        # Commit the session to persist changes to the database
        db.session.commit()

        # get updated sample's variants
        sample: Samples = db.session.query(Samples).get(sample_id)
        variants_samples: List[VariantsSamples] = sample.variants_samples
        updated_variants, updated_not_sample_variants = get_sample_variants(variants_samples)

        return InternalResponse({'isSuccess': True, 'isSampleDeleted': False, 'updatedVariants': updated_variants, "updatedNotSampleVariants": updated_not_sample_variants}, 200)
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since addition of new variants to sample {sample_id} in DB failed due to error: {e}')
        return InternalResponse({'isSuccess': False}, 500)


def add_new_sample_to_db(sample_id: str, phenotypes: List) -> Samples:
    # check if sample already exists
    sample: Samples = db.session.query(Samples).filter(Samples.id == sample_id).one_or_none()

    # if the sample is new, add it to database
    if sample is None:
        sample = Samples(id=sample_id, genome_version='GRCh37')
        db.session.add(sample)

        db.session.flush()

    sample_ontology_term_ids = [o.ontology_term_id for o in sample.ontology_term]

    # append phenotypes to the respective sample, if the sample does not already have that phenotype
    for phenotype_term in phenotypes:
        if phenotype_term['ontologyId'] not in sample_ontology_term_ids:
            append_phenotype_to_sample(sample, phenotype_term)

    return sample
