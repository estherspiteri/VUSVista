from datetime import datetime
from typing import List, Dict

import pandas as pd
from sqlalchemy import desc

from server import db
from server.helpers.data_helper import convert_df_to_list
from server.models import ExternalReferences, Variants, DbSnp, Clinvar, VariantsSamples, Genotype, Samples
from server.services.clinvar_service import get_last_saved_clinvar_update


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


def retrieve_vus_from_db(vus_id: int) -> Dict:
    variant: Variants = db.session.query(Variants).filter(Variants.id == vus_id).first()

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

    # retrieve samples that have this variant
    variant_samples: List[VariantsSamples] = db.session.query(VariantsSamples).filter(VariantsSamples.variant_id == variant.id).all()
    variant_data['samples'] = [{'id': vs.sample_id, 'hgvs': vs.variant_hgvs.hgvs} for vs in variant_samples]

    # retrieve all the unique phenotypes that these samples have
    phenotypes = []
    phenotype_ids = []

    samples: List[Samples] = [vs.sample for vs in variant_samples]
    for s in samples:
        for term in s.ontology_term:
            if term.ontology_term_id not in phenotype_ids:
                phenotype_ids.append(term.ontology_term_id)
                phenotypes.append({'ontologyId': term.ontology_term_id, 'name': term.term_name})

    variant_data['phenotypes'] = phenotypes

    return variant_data
