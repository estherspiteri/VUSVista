from datetime import datetime
from typing import List

import pandas as pd

from server import db
from server.helpers.data_helper import convert_df_to_list
from server.models import ExternalReferences, Variants, DbSnp, Clinvar, VariantsSamples, Genotype


# TODO: merge with the function found in vus_preprocess_service.py
def add_missing_columns(vus_df: pd.DataFrame) -> pd.DataFrame:
    # insert columns for clinvar clinical significance and error messages
    vus_df['clinvarClassification'] = ""
    vus_df['clinvarClassificationLastEval'] = ""
    vus_df['clinvarClassificationReviewStatus'] = ""
    vus_df['clinvarErrorMsg'] = ""
    vus_df['clinvarCanonicalSpdi'] = ""
    vus_df['clinvarUid'] = ""

    # insert columns for dbsnp and error messages
    vus_df['rsid'] = ""
    vus_df['rsidDbsnpVerified'] = False
    vus_df['rsidDbsnpErrorMsgs'] = ""

    # create new column to determine whether the variant has already been stored in the db or not
    vus_df.insert(0, 'Exists in DB', False)

    return vus_df


def retrieve_all_vus_from_db():
    variants: List[Variants] = db.session.query(Variants).all()

    variants_data = [{'variantId': v.variant_id, 'chromosome': v.chromosome,
                      'chromosomePosition': v.chromosome_position, 'gene': v.gene_name,
                      'type': v.variant_type.value, 'refAllele': v.ref, 'observedAllele': v.alt,
                      'classification': v.classification.value} for v in variants]

    # store the variants into a dataframe
    vus_df = pd.DataFrame(variants_data)

    vus_df = add_missing_columns(vus_df)

    vus_df_copy = vus_df.copy()

    # iterate through the dataframe
    for index, row in vus_df_copy.iterrows():
        # retrieve all external references related to that variant
        external_references: List[ExternalReferences] = db.session.query(ExternalReferences).filter(
            ExternalReferences.variant_id == row['variantId']
        ).all()

        for ref in external_references:
            if ref.db_type == 'db_snp':
                # retrieve dbsnp entry related to the variant
                dbsnp: DbSnp = db.session.query(DbSnp).filter(
                    DbSnp.external_references_id == ref.external_references_id
                ).one_or_none()

                vus_df.at[index, 'rsid'] = dbsnp.db_snp_id
                vus_df.at[index, 'rsidDbsnpVerified'] = len(ref.error_msg) == 0
                vus_df.at[index, 'rsidDbsnpErrorMsgs'] = ref.error_msg

            elif ref.db_type == 'clinvar':
                # retrieve clinvar entry related to the variant
                clinvar: Clinvar = db.session.query(Clinvar).filter(
                    Clinvar.external_references_id == ref.external_references_id
                ).one_or_none()

                if clinvar.last_evaluated is not None:
                    clinvar_last_evaluated = datetime.strftime(clinvar.last_evaluated, '%Y/%m/%d %H:%M')
                else:
                    clinvar_last_evaluated = None

                # populate the clinvar fields
                vus_df.at[index, 'clinvarUid'] = clinvar.clinvar_id
                vus_df.at[index, 'clinvarCanonicalSpdi'] = clinvar.canonical_spdi
                vus_df.at[index, 'clinvarClassification'] = clinvar.classification
                vus_df.at[index, 'clinvarClassificationReviewStatus'] = clinvar.review_status
                vus_df.at[index, 'clinvarClassificationLastEval'] = clinvar_last_evaluated
                vus_df.at[index, 'clinvarErrorMsg'] = ref.error_msg

        # retrieve all samples related to that variant
        variant_samples: List[VariantsSamples] = (db.session.query(VariantsSamples)
                                                  .filter(VariantsSamples.variant_id == row['variantId'])).all()

        num_heterozygous = len([s for s in variant_samples if s.genotype == Genotype.HETEROZYGOUS])
        num_homozygous = len(variant_samples) - num_heterozygous

        vus_df.at[index, 'numHeterozygous'] = num_heterozygous
        vus_df.at[index, 'numHomozygous'] = num_homozygous

    var_list = convert_df_to_list(vus_df)

    return var_list
