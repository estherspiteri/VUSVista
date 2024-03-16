from typing import List, Dict

from flask import current_app
from sqlalchemy.exc import SQLAlchemyError

from server import db
from server.models import AcmgRules, VariantsSamplesAcmgRules


def get_acmg_rules() -> List[Dict[str, str]]:
    acmg_rules: List[AcmgRules] = db.session.query(AcmgRules).all()

    return [{'id': r.id, 'name': r.rule_name.value} for r in acmg_rules]


def add_acmg_rule_to_sample_variant(sample_id: str, variant_id: int, acmg_rule_id: str):
    acmg_rule: AcmgRules = db.session.query(AcmgRules).filter(AcmgRules.id == acmg_rule_id).first()
    variants_samples_acmg_rules = VariantsSamplesAcmgRules(sample_id=sample_id, variant_id=variant_id, acmg_rule_id=acmg_rule_id, rule_name=acmg_rule.rule_name)

    db.session.add(variants_samples_acmg_rules)

    try:
        # Commit the session to persist changes to the database
        db.session.commit()
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since insertion of SamplesVariantsAcmgRules entry in DB failed due to error: {e}')


def remove_acmg_rule_to_sample_variant(sample_id: str, variant_id: int, acmg_rule_id: str):
    samples_variants_acmg_rules = db.session.query(VariantsSamplesAcmgRules).filter_by(sample_id=sample_id, variant_id=variant_id, acmg_rule_id=acmg_rule_id).first()

    db.session.delete(samples_variants_acmg_rules)

    try:
        # Commit the session to persist changes to the database
        db.session.commit()
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since deletion of SamplesVariantsAcmgRules entry in DB failed due to error: {e}')
