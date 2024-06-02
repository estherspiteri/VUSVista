from typing import List, Dict

from flask import current_app
from sqlalchemy.exc import SQLAlchemyError

from server import db
from server.models import AcmgRules, VariantsAcmgRules
from server.responses.internal_response import InternalResponse


def get_acmg_rules() -> List[Dict[str, str]]:
    acmg_rules: List[AcmgRules] = db.session.query(AcmgRules).all()

    return [{'id': r.id, 'name': r.rule_name.value, 'description': r.description, 'defaultStrength': r.default_strength.value} for r in acmg_rules]


def add_acmg_rule_to_variant(variant_id: int, acmg_rule_id: str):
    acmg_rule: AcmgRules = db.session.query(AcmgRules).filter(AcmgRules.id == acmg_rule_id).first()
    variants_acmg_rules = VariantsAcmgRules(variant_id=variant_id, acmg_rule_id=acmg_rule_id, rule_name=acmg_rule.rule_name)

    db.session.add(variants_acmg_rules)


def remove_acmg_rule_from_variant(variant_id: int, acmg_rule_id: str):
    variants_acmg_rules = db.session.query(VariantsAcmgRules).filter_by(variant_id=variant_id,
                                                                        acmg_rule_id=acmg_rule_id).first()

    db.session.delete(variants_acmg_rules)

