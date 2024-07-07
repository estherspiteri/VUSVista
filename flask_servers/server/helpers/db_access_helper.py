import html
from typing import List

from flask import current_app
from pandas import Series
from sqlalchemy.exc import MultipleResultsFound

from server import db
from server.models import Variants


def get_variant_from_db(vus_df_row: Series) -> Variants | None:
    try:
        # retrieve variant if it already exists in DB
        variant: Variants = db.session.query(Variants).filter(
            Variants.chromosome == vus_df_row['Chr'],
            Variants.chromosome_position == vus_df_row['Position'],
            Variants.variant_type == html.unescape(vus_df_row['Type']),
            Variants.ref == vus_df_row['Reference'],
            Variants.alt == vus_df_row['Alt'],
            Variants.gene_id == vus_df_row['Gene Id']
        ).one_or_none()

        return variant
    except MultipleResultsFound as e:
        current_app.logger.error(
            f'The following variant :\n{vus_df_row}\n exists more than once in db!', e)
        return None
