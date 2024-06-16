from typing import List, Dict, Tuple

import pandas as pd
from flask import current_app
from sqlalchemy.exc import SQLAlchemyError

from server import db
from server.helpers.data_helper import convert_df_to_list
from server.models import ExternalReferences, Variants, DbSnp, Clinvar, VariantsSamples, Genotype, Samples, Phenotypes, \
    FileUploads, Publications, AutoClinvarUpdates, VariantHgvs
from server.responses.internal_response import InternalResponse
from server.services.clinvar_service import get_last_saved_clinvar_update
from server.services.consequence_service import get_consequences_for_new_vus
from server.services.phenotype_service import append_phenotype_to_sample
from server.services.samples_service import add_new_sample_to_db
from server.services.variants_samples_service import store_upload_details_for_variant_sample, add_variant_sample_to_db


def retrieve_all_vus_summaries_from_db():
    variants: List[Variants] = db.session.query(Variants).all()

    # sort by id
    variants.sort(key=lambda x: x.id)

    variants_data = [{'id': v.id, 'chromosome': v.chromosome,
                      'chromosomePosition': v.chromosome_position, 'gene': v.gene_name,
                      'refAllele': v.ref, 'altAllele': v.alt} for v in variants]

    # store the variants into a dataframe
    vus_df = pd.DataFrame(variants_data)

    # insert columns for dbsnp
    vus_df['rsid'] = ""
    vus_df['rsidDbsnpVerified'] = False

    vus_df_copy = vus_df.copy()

    # iterate through the dataframe
    for index, row in vus_df_copy.iterrows():
        # retrieve all external references related to that variant
        external_references: List[ExternalReferences] = db.session.query(ExternalReferences).filter(
            ExternalReferences.variant_id == row['id']
        ).all()

        for ref in external_references:
            if ref.db_type == 'db_snp':
                # retrieve dbsnp entry related to the variant
                dbsnp: DbSnp = db.session.query(DbSnp).filter(
                    DbSnp.external_db_snp_id == ref.id
                ).one_or_none()

                vus_df.at[index, 'rsid'] = dbsnp.rsid
                vus_df.at[index, 'rsidDbsnpVerified'] = len(ref.error_msg) == 0

    var_list = convert_df_to_list(vus_df)

    return var_list


def get_variant_samples(variant_samples: List[VariantsSamples]) -> Tuple[List, List]:
    variant_samples_list = [
        {'id': vs.sample_id, 'hgvs': vs.variant_hgvs.hgvs, 'noOfVariants': len(vs.sample.variants_samples), 'consequence': vs.consequence} for vs in
        variant_samples]

    # retrieve samples that do not have this variant
    samples_ids = [vs.sample_id for vs in variant_samples]
    not_variants_samples: List[Samples] = db.session.query(Samples).filter(Samples.id.not_in(samples_ids)).all()
    not_variant_samples_list = [{'id': s.id, 'noOfVariants': len(s.variants_samples)} for s in
                                     not_variants_samples]

    return variant_samples_list, not_variant_samples_list


def get_variant_phenotypes_from_db(samples: List[Samples]):
    # retrieve all the unique phenotypes that these samples have
    phenotypes = []
    phenotype_ids = []

    for s in samples:
        for term in s.ontology_term:
            if term.ontology_term_id not in phenotype_ids:
                phenotype_ids.append(term.ontology_term_id)
                phenotypes.append({'ontologyId': term.ontology_term_id, 'name': term.term_name})

    return phenotypes


def retrieve_vus_from_db(vus_id: int) -> (Dict | None):
    variant: Variants = db.session.query(Variants).filter(Variants.id == vus_id).one_or_none()

    if variant is None:
        return None

    variant_data = {'id': variant.id, 'chromosome': variant.chromosome,
                    'chromosomePosition': variant.chromosome_position, 'gene': variant.gene_name,
                    'type': variant.variant_type.value, 'refAllele': variant.ref, 'altAllele': variant.alt,
                    'classification': variant.classification.value,
                    'acmgRuleIds': [r.acmg_rule_id for r in variant.variants_acmg_rules], 'numOfPublications': len(variant.variants_publications)}

    # retrieve all external references related to that variant
    external_references: List[ExternalReferences] = db.session.query(ExternalReferences).filter(
        ExternalReferences.variant_id == variant.id
    ).all()

    for ref in external_references:
        if ref.db_type == 'db_snp':
            # retrieve dbsnp entry related to the variant
            dbsnp: DbSnp = db.session.query(DbSnp).filter(
                DbSnp.external_db_snp_id == ref.id
            ).one_or_none()

            variant_data['rsid'] = dbsnp.rsid
            variant_data['rsidDbsnpVerified'] = len(ref.error_msg) == 0
            variant_data['rsidDbsnpErrorMsgs'] = ref.error_msg

        elif ref.db_type == 'clinvar':
            # retrieve clinvar entry related to the variant
            clinvar: Clinvar = db.session.query(Clinvar).filter(
                Clinvar.external_clinvar_id == ref.id
            ).one_or_none()

            if clinvar is not None:
                auto_clinvar_update_id, review_status, classification, last_evaluated = get_last_saved_clinvar_update(clinvar.id)

                # populate the clinvar fields
                variant_data['clinvarId'] = clinvar.id
                variant_data['clinvarVariationId'] = clinvar.variation_id
                variant_data['clinvarCanonicalSpdi'] = clinvar.canonical_spdi
                variant_data['clinvarClassification'] = classification
                variant_data['clinvarClassificationReviewStatus'] = review_status
                variant_data['clinvarClassificationLastEval'] = last_evaluated
                variant_data['clinvarErrorMsg'] = ref.error_msg

    # retrieve all samples related to that variant
    variant_samples: List[VariantsSamples] = (db.session.query(VariantsSamples)
                                              .filter(VariantsSamples.variant_id == variant.id)).all()

    num_heterozygous = len([s for s in variant_samples if s.genotype == Genotype.HETEROZYGOUS])
    num_homozygous = len(variant_samples) - num_heterozygous

    variant_data['numHeterozygous'] = num_heterozygous
    variant_data['numHomozygous'] = num_homozygous

    variant_data['samples'], variant_data['notVusSamples'] = get_variant_samples(variant_samples)

    samples: List[Samples] = [vs.sample for vs in variant_samples]

    # retrieve all the unique phenotypes that these samples have
    variant_data['phenotypes'] = get_variant_phenotypes_from_db(samples)

    return variant_data


def delete_variant_entry(variant_id: str) -> InternalResponse:
    variant: Variants = db.session.query(Variants).get(variant_id)

    if not variant:
        current_app.logger.error(f'Variant with ID: {variant_id} not found')
        return InternalResponse({'isSuccess': False}, 404)

    current_app.logger.info(
        f"Deleting variant with ID: {variant_id}")
    db.session.delete(variant)

    # deleting any publications without variants_publications
    publications_query = db.session.query(Publications).filter(~Publications.variants_publications.any())
    publications_query.delete()

    # deleting any auto_clinvar_updates without auto_clinvar_eval_dates
    auto_clinvar_updates_query = db.session.query(AutoClinvarUpdates).filter(AutoClinvarUpdates.auto_clinvar_eval_dates == None)
    auto_clinvar_updates_query.delete()

    # deleting any samples without variants_samples
    samples_query = db.session.query(Samples).filter(~Samples.variants_samples.any())
    samples_query.delete()

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
            f'Rollback carried out since deletion of variant {variant_id} in DB failed due to error: {e}')
        return InternalResponse({'isSuccess': False}, 500)


def commit_samples_update_to_variant(variant_id: int):
    try:
        # Commit the session to persist changes to the database
        db.session.commit()

        # retrieve all samples related to that variant
        variant_samples: List[VariantsSamples] = (db.session.query(VariantsSamples)
                                                  .filter(VariantsSamples.variant_id == variant_id)).all()

        # get updated variant's samples
        updated_samples, updated_not_variant_samples = get_variant_samples(variant_samples)

        samples: List[Samples] = [vs.sample for vs in variant_samples]
        updated_phenotypes = get_variant_phenotypes_from_db(samples)

        return InternalResponse({'isSuccess': True, 'updatedSamples': updated_samples, "updatedNotVariantSamples": updated_not_variant_samples, "updatedPhenotypes": updated_phenotypes}, 200)
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since samples update for variant {variant_id} in DB failed due to error: {e}')
        return InternalResponse({'isSuccess': False}, 500)


def add_samples_to_variant(variant_id: int, samples_to_add: List) -> InternalResponse:
    for s in samples_to_add:
        hgvs: VariantHgvs = db.session.query(VariantHgvs).filter(VariantHgvs.variant_id == variant_id, VariantHgvs.hgvs == s['hgvs']).one_or_none()

        if hgvs is None:
            hgvs = VariantHgvs(variant_id=variant_id, hgvs=s['hgvs'], is_updated=False)
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

        variant_sample = VariantsSamples(variant_id=variant_id, sample_id=s['sampleId'], genotype=s['genotype'].upper(), variant_hgvs_id=hgvs.id, consequence=hgvs_consequence_dict.get(hgvs.hgvs, ""))
        db.session.add(variant_sample)

        store_upload_details_for_variant_sample(None, False, s['sampleId'], variant_id)

        if 'phenotypes' in s.keys():
            sample: Samples = db.session.query(Samples).filter(Samples.id == s['sampleId']).first()
            sample_ontology_term_ids = [o.ontology_term_id for o in sample.ontology_term]
            for p in s['phenotypes']:
                if p['ontologyId'] not in sample_ontology_term_ids:
                    append_phenotype_to_sample(sample, p)

    return commit_samples_update_to_variant(variant_id)


def add_new_sample_to_variant(variant_id: int, sample_to_add: Dict) -> InternalResponse:
    phenotypes = []
    if "phenotypes" in sample_to_add.keys():
        phenotypes = sample_to_add["phenotypes"]

    new_sample = add_new_sample_to_db(sample_to_add["sampleId"], phenotypes)

    db.session.flush()

    hgvs_consequence_dict = {}

    # get variant consequences though HGVS
    get_consequences_for_new_vus_res = get_consequences_for_new_vus([sample_to_add["hgvs"]])

    if get_consequences_for_new_vus_res.status != 200:
        current_app.logger.error(
            f'Retrieval of variant consequence failed 500')
    else:
        hgvs_consequence_dict = get_consequences_for_new_vus_res.data['consequences_dict']

    # store the variants sample
    add_variant_sample_to_db(variant_id, new_sample.id, sample_to_add["hgvs"], sample_to_add["genotype"], hgvs_consequence_dict.get(sample_to_add["hgvs"], ""))

    # store the upload details related to this variant & sample
    store_upload_details_for_variant_sample(None, False, new_sample.id, variant_id)

    return commit_samples_update_to_variant(variant_id)


def remove_sample_from_variant(variant_id: int, sample_ids_to_remove: List[str]):
    db.session.query(VariantsSamples).filter(VariantsSamples.variant_id == variant_id,
                                             VariantsSamples.sample_id.in_(
                                                 sample_ids_to_remove)).delete()

    # deleting any samples without variants_samples
    db.session.query(Samples).filter(~Samples.variants_samples.any()).delete()

    # deleting any file_uploads without variants_samples_uploads
    db.session.query(FileUploads).filter(~FileUploads.variants_samples_uploads.any()).delete()

    # deleting any phenotypes without samples
    db.session.query(Phenotypes).filter(~Phenotypes.sample.any()).delete()

    return commit_samples_update_to_variant(variant_id)
