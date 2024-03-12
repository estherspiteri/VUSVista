from typing import List

from flask import current_app
from sqlalchemy.exc import SQLAlchemyError

from server import db
from server.models import AcmgRules, SamplesVariantsAcmgRules


def get_acmg_rule_names() -> List[str]:
    acmg_rules: List[AcmgRules] = db.session.query(AcmgRules.rule_name).all()

    return [r.rule_name.value for r in acmg_rules]


def add_acmg_rule_to_sample_variant(sample_id: str, variant_id: int, rule_name: str):
    samples_variants_acmg_rules = SamplesVariantsAcmgRules(sample_id=sample_id, variant_id=variant_id, rule_name=rule_name)

    db.session.add(samples_variants_acmg_rules)

    try:
        # Commit the session to persist changes to the database
        db.session.commit()
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since insertion of SamplesVariantsAcmgRules entry in DB failed due to error: {e}')


def remove_acmg_rule_to_sample_variant(sample_id: str, variant_id: int, rule_name: str):
    samples_variants_acmg_rules = db.session.query(SamplesVariantsAcmgRules).filter_by(sample_id=sample_id, variant_id=variant_id, rule_name=rule_name).first()

    db.session.delete(samples_variants_acmg_rules)

    try:
        # Commit the session to persist changes to the database
        db.session.commit()
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since deletion of SamplesVariantsAcmgRules entry in DB failed due to error: {e}')
