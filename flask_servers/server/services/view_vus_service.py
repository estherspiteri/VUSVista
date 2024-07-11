from typing import List, Dict, Tuple

import pandas as pd
from flask import current_app
from flask_login import current_user
from sqlalchemy import desc
from sqlalchemy.exc import SQLAlchemyError

from server import db
from server.helpers.data_helper import convert_df_to_list
from server.models import ExternalReferences, Variants, DbSnp, Clinvar, VariantsSamples, Genotype, Samples, Phenotypes, \
    FileUploads, Publications, AutoClinvarUpdates, VariantHgvs, VariantsPublications, AutoPublicationEvalDates, \
    AutoClinvarEvalDates, VariantsSamplesUploads
from server.responses.internal_response import InternalResponse
from server.services.clinvar_service import get_last_saved_clinvar_update, store_clinvar_info, clinvar_clinical_significance_pipeline_single
from server.services.consequence_service import get_consequences_for_new_vus
from server.services.phenotype_service import append_phenotype_to_sample
from server.services.publications_service import update_variant_publications
from server.services.samples_service import add_new_sample_to_db
from server.services.variants_samples_service import store_upload_details_for_variant_sample, add_variant_sample_to_db


def retrieve_vus_summaries_from_db(variants: List[Variants]):
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
    vus_df['rsidReviewRequired'] = False

    # insert column for clinvar flag
    vus_df['isFoundInClinvar'] = False

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

                # only return RSID if it is valid for the respective variant
                if len(ref.error_msg) == 0:
                    vus_df.at[index, 'rsid'] = dbsnp.rsid
                else:
                    vus_df.at[index, 'rsidReviewRequired'] = True

            # only return True Clinvar flag if it is valid for the respective variant
            elif ref.db_type == 'clinvar' and len(ref.error_msg) == 0:
                vus_df.at[index, 'isFoundInClinvar'] = len(ref.error_msg) == 0

    var_list = convert_df_to_list(vus_df)

    return var_list


def get_variant_samples(variant_samples: List[VariantsSamples]) -> Tuple[List, List]:
    variant_samples_list = []
    for vs in variant_samples:
        hgvs = None
        if vs.variant_hgvs is not None:
            hgvs = vs.variant_hgvs.hgvs

        variant_samples_list.append({'id': vs.sample_id, 'hgvs': hgvs, 'noOfVariants': len(vs.sample.variants_samples),
                                     'consequence': vs.consequence})

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
                    'type': variant.variant_type, 'refAllele': variant.ref, 'altAllele': variant.alt,
                    'classification': variant.classification.value,
                    'acmgRuleIds': [r.acmg_rule_id for r in variant.variants_acmg_rules],
                    'numOfPublications': len(variant.variants_publications)}

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
                auto_clinvar_update_id, review_status, classification, last_evaluated = get_last_saved_clinvar_update(
                    clinvar.id)

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
    auto_clinvar_updates_query = db.session.query(AutoClinvarUpdates).filter(
        AutoClinvarUpdates.auto_clinvar_eval_dates == None)
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

        return InternalResponse({'isSuccess': True, 'updatedSamples': updated_samples,
                                 "updatedNotVariantSamples": updated_not_variant_samples,
                                 "updatedPhenotypes": updated_phenotypes}, 200)
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since samples update for variant {variant_id} in DB failed due to error: {e}')
        return InternalResponse({'isSuccess': False}, 500)


def add_samples_to_variant(variant_id: int, samples_to_add: List) -> InternalResponse:
    for s in samples_to_add:
        hgvs_id = None
        consequence = None

        if 'hgvs' in s.keys():
            hgvs: VariantHgvs = db.session.query(VariantHgvs).filter(VariantHgvs.variant_id == variant_id,
                                                                     VariantHgvs.hgvs == s['hgvs']).one_or_none()

            if hgvs is None:
                hgvs = VariantHgvs(variant_id=variant_id, hgvs=s['hgvs'], is_updated=False)
                db.session.add(hgvs)
                db.session.flush()

            hgvs_id = hgvs.id

            # get variant consequences though HGVS
            get_consequences_for_new_vus_res = get_consequences_for_new_vus([hgvs.hgvs])

            if get_consequences_for_new_vus_res.status != 200:
                current_app.logger.error(
                    f'Retrieval of variant consequence failed 500')
            else:
                hgvs_consequence_dict = get_consequences_for_new_vus_res.data['consequences_dict']
                consequence = hgvs_consequence_dict.get(hgvs.hgvs, "")

        variant_sample = VariantsSamples(variant_id=variant_id, sample_id=s['sampleId'], genotype=s['genotype'].upper(),
                                         variant_hgvs_id=hgvs_id, consequence=consequence)
        db.session.add(variant_sample)

        store_upload_details_for_variant_sample(None, False, s['sampleId'], variant_id, current_user.id)

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

    hgvs = None
    consequence = None
    if 'hgvs' in sample_to_add.keys():
        hgvs = sample_to_add["hgvs"]
        # get variant consequences though HGVS
        get_consequences_for_new_vus_res = get_consequences_for_new_vus([hgvs])

        if get_consequences_for_new_vus_res.status != 200:
            current_app.logger.error(
                f'Retrieval of variant consequence failed 500')
        else:
            hgvs_consequence_dict = get_consequences_for_new_vus_res.data['consequences_dict']
            consequence = hgvs_consequence_dict.get(hgvs, "")

    # store the variants sample
    add_variant_sample_to_db(variant_id, new_sample.id, hgvs, sample_to_add["genotype"],
                             consequence)

    # store the upload details related to this variant & sample
    store_upload_details_for_variant_sample(None, False, new_sample.id, variant_id, current_user.id)

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


def update_variant_rsid(variant_id: int, new_rsid: str):
    variant: Variants = db.session.query(Variants).get(variant_id)

    db_snp: DbSnp | None = None
    db_snp_ref: ExternalReferences | None = None

    clinvar: Clinvar | None = None
    clinvar_ref: ExternalReferences | None = None
    clinvar_error_msg = ''
    clinvar_clinical_significance = None

    # check if DbSNP entry exists for given variant
    for ref in variant.external_references:
        if ref.db_type == 'db_snp':
            ref.db_snp.rsid = new_rsid
            db_snp = ref.db_snp
            db_snp_ref = ref
            db_snp_ref.error_msg = ""
        else:
            clinvar = ref.clinvar
            clinvar_ref = ref

    # create new DbSNP entry
    if db_snp is None:
        db_snp_ref = ExternalReferences(variant_id=variant_id,
                                        db_type='db_snp',
                                        error_msg=None)
        db.session.add(db_snp_ref)
        db.session.flush()

        db_snp = DbSnp(rsid=new_rsid,
                       external_db_snp_id=db_snp_ref.id)
        db.session.add(db_snp)

    # remove all automatically added publications since they are related to old RSID
    db.session.query(VariantsPublications).filter(VariantsPublications.variant_id == variant_id, VariantsPublications.is_manually_added.is_(False)).delete()

    manual_publications = [p for p in variant.variants_publications if p.is_manually_added]
    variant.variants_publications = manual_publications

    db.session.flush()

    # remove publications which are no longer relate dto a variant
    db.session.query(Publications).filter(~Publications.variants_publications.any()).delete()

    # remove publication eval dates for variant
    db.session.query(AutoPublicationEvalDates).filter(AutoPublicationEvalDates.variant_id == variant_id).delete()

    db.session.flush()

    # get publications related to new RSID
    hgvs = None
    # get one of the variant's HGVS
    if len(variant.variant_hgvs) > 0:
        hgvs = variant.variant_hgvs[0].hgvs.split(' ')[0]

    update_variant_publications(variant, hgvs, new_rsid)

    clinvar_clinical_significance_pipeline_res = clinvar_clinical_significance_pipeline_single('GRCh37',
                                                                                        new_rsid,
                                                                                        variant.gene_name,
                                                                                        variant.chromosome,
                                                                                        variant.chromosome_position)

    if clinvar_clinical_significance_pipeline_res.status != 200:
        current_app.logger.error(
            f"ClinVar clinical significance pipeline failed for variant with RSID {new_rsid}!")
    else:
        # execute pipeline
        is_success, clinical_significance, canonical_spdi, variation_id, error_msg = (
            clinvar_clinical_significance_pipeline_res.data)

        clinvar_error_msg = error_msg
        clinvar_clinical_significance = clinical_significance

        # if clinvar entry found
        if len(variation_id) > 0:
            if clinvar is None:
                clinvar_ref = ExternalReferences(variant_id=variant_id,
                                                 db_type='clinvar',
                                                 error_msg=error_msg)
                db.session.add(clinvar_ref)
                db.session.flush()

                clinvar = Clinvar(variation_id=variation_id,
                                  external_clinvar_id=clinvar_ref.id,
                                  canonical_spdi=canonical_spdi)
                db.session.add(clinvar)
                db.session.flush()
            else:
                # delete existing entries of auto_eval_dates and auto_clinvar_updates
                auto_clinvar_eval_dates_query = db.session.query(AutoClinvarEvalDates).filter(
                    AutoClinvarEvalDates.clinvar_id == clinvar.id)
                auto_clinvar_eval_dates: List[AutoClinvarEvalDates] = auto_clinvar_eval_dates_query.all()
                for eval_dates in auto_clinvar_eval_dates:
                    if eval_dates.auto_clinvar_update_id is not None:
                        db.session.delete(eval_dates.auto_clinvar_update)

                auto_clinvar_eval_dates_query.delete()

                clinvar.variation_id = variation_id
                clinvar.canonical_spdi = canonical_spdi
                clinvar_ref.error_msg = error_msg

            store_clinvar_info(clinvar.id, clinical_significance['description'],
                               clinical_significance['review_status'], clinical_significance['last_evaluated'],
                               True)
        # if no clinvar entry found and the variant used to have a clinvar entry
        elif clinvar is not None:
            # delete existing entries of auto_eval_dates and auto_clinvar_updates
            auto_clinvar_eval_dates_query = db.session.query(AutoClinvarEvalDates).filter(
                AutoClinvarEvalDates.clinvar_id == clinvar.id)
            auto_clinvar_eval_dates: List[AutoClinvarEvalDates] = auto_clinvar_eval_dates_query.all()
            for eval_dates in auto_clinvar_eval_dates:
                if eval_dates.auto_clinvar_update_id is not None:
                    db.session.delete(eval_dates.auto_clinvar_update)

            auto_clinvar_eval_dates_query.delete()

            db.session.delete(clinvar)

            db.session.delete(clinvar_ref)

            clinvar = None
    try:
        # Commit the session to persist changes to the database
        db.session.commit()

        updated_external_ref_data = {'rsid': new_rsid, 'rsidDbsnpVerified': len(db_snp_ref.error_msg) == 0,
                                     'rsidDbsnpErrorMsgs': db_snp_ref.error_msg, 'numOfPublications': len(variant.variants_publications)}

        if clinvar is not None:
            # populate the clinvar fields
            updated_external_ref_data['clinvarId'] = clinvar.id
            updated_external_ref_data['clinvarVariationId'] = clinvar.variation_id
            updated_external_ref_data['clinvarCanonicalSpdi'] = clinvar.canonical_spdi
            if len(clinvar_clinical_significance.keys()) > 0:
                updated_external_ref_data['clinvarClassification'] = clinvar_clinical_significance['description']
                updated_external_ref_data['clinvarClassificationReviewStatus'] = clinvar_clinical_significance[
                    'review_status']
                updated_external_ref_data['clinvarClassificationLastEval'] = clinvar_clinical_significance['last_evaluated']
            else:
                updated_external_ref_data['clinvarClassification'] = ""
                updated_external_ref_data['clinvarClassificationReviewStatus'] = ""
                updated_external_ref_data['clinvarClassificationLastEval'] = ""
            updated_external_ref_data['clinvarErrorMsg'] = clinvar_error_msg

        return InternalResponse({'isSuccess': True, 'updated_external_ref_data': updated_external_ref_data}, 200)
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since RSID update for variant ID {variant_id} to {new_rsid} in DB failed due to error: {e}')
        return InternalResponse({'isSuccess': False}, 500)


def get_latest_added_vus(num_of_variants: int):
    latest_variant_sample_uploads: List[VariantsSamplesUploads] = db.session.query(VariantsSamplesUploads).order_by(
        desc(VariantsSamplesUploads.date_uploaded)).all()

    updated_num_variants = num_of_variants

    latest_variant_ids = list(set([u.variant_id for u in latest_variant_sample_uploads]))
    if len(latest_variant_ids) < num_of_variants:
        updated_num_variants = len(latest_variant_ids)

    latest_variant_ids = latest_variant_ids[:updated_num_variants]

    variants: List[Variants] = db.session.query(Variants).filter(Variants.id.in_(latest_variant_ids)).all()

    var_list = retrieve_vus_summaries_from_db(variants)

    return InternalResponse({'var_list': var_list}, 500)
