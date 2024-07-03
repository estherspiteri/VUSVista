import math

import pandas as pd
from typing import List, Dict

from sqlalchemy.orm import DeclarativeMeta

from server import db
from server.models import Variants, ExternalReferences, DbSnp


def alchemy_encoder(obj):
    """JSON encoder function for SQLAlchemy special classes."""
    if isinstance(obj.__class__, DeclarativeMeta):
        # an SQLAlchemy class
        # this will give the individual attributes, in case you want to exclude some fields you can do so here
        data = {}
        for column in obj.__table__.columns:
            data[column.name] = getattr(obj, column.name)
        # for relation in obj.__mapper__.relationships:
        #     data[relation.key] = [alchemy_encoder(i) for i in getattr(obj, relation.key)]
        return data
    raise TypeError(f"Object of type '{obj.__class__.__name__}' is not JSON serializable")


def convert_df_to_list(df: pd.DataFrame) -> List:
    # convert dataframe to list
    final_list = []
    for index, row in df.iterrows():
        final_list.append(row.to_dict())
    return final_list


def prep_unprocessed_vus_dict_for_react(vus: Dict) -> Dict:
    ref = vus['Reference']
    if isinstance(ref, float) and math.isnan(ref):
        ref = None

    alt = vus['Alt']
    if isinstance(alt, float) and math.isnan(alt):
        alt = None

    # to match React camelCase syntax
    new_vus = {'locus': vus['Locus'], 'type': vus['Type'], 'genotype': vus['Genotype'], 'refAllele': ref,
               'altAllele': alt}

    return new_vus


def prep_vus_df_for_react(vus_df: pd.DataFrame) -> pd.DataFrame:
    # to match React camelCase syntax
    new_vus_df = pd.DataFrame()
    new_vus_df['id'] = vus_df['Variant Id']
    new_vus_df['chromosome'] = vus_df['Chr']
    new_vus_df['chromosomePosition'] = vus_df['Position']
    new_vus_df['gene'] = vus_df['Gene']
    new_vus_df['type'] = vus_df['Type']
    new_vus_df['genotype'] = vus_df['Genotype']
    new_vus_df['refAllele'] = vus_df['Reference']
    new_vus_df['altAllele'] = vus_df['Alt']
    new_vus_df['classification'] = vus_df['Classification']
    new_vus_df['rsid'] = vus_df['RSID']
    new_vus_df['rsidDbsnpVerified'] = vus_df['RSID dbSNP verified']
    new_vus_df['rsidDbsnpErrorMsgs'] = vus_df['RSID dbSNP errorMsgs']
    new_vus_df['clinvarErrorMsg'] = vus_df['Clinvar error msg']
    new_vus_df['clinvarClassification'] = vus_df['Clinvar classification']
    new_vus_df['clinvarClassificationLastEval'] = vus_df['Clinvar classification last eval']
    new_vus_df['clinvarClassificationReviewStatus'] = vus_df['Clinvar classification review status']
    new_vus_df['clinvarCanonicalSpdi'] = vus_df['Clinvar canonical spdi']
    new_vus_df['clinvarVariationId'] = vus_df['Clinvar variation id']

    return new_vus_df


def get_variant_summary(variant: Variants, include_rsid=False) -> Dict:
    variant_summary = {'id': variant.id, 'chromosome': variant.chromosome, 'chromosomePosition': variant.chromosome_position,
                       'gene': variant.gene_name, 'altAllele': variant.alt, 'refAllele': variant.ref}

    if include_rsid:
        db_snp_external_ref: ExternalReferences = db.session.query(ExternalReferences).filter(
            ExternalReferences.variant_id == variant.id, ExternalReferences.db_type == 'db_snp'
        ).one_or_none()

        if db_snp_external_ref is not None:
            db_snp: DbSnp = db.session.query(DbSnp).filter(
                DbSnp.external_db_snp_id == db_snp_external_ref.id
            ).first()

            variant_summary['rsid'] = db_snp.rsid

    return variant_summary
