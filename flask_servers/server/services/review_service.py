from datetime import datetime
from typing import List, Dict, Tuple

from flask import current_app
from flask_login import current_user
from sqlalchemy.exc import SQLAlchemyError

from server import db
from server.helpers.data_helper import get_variant_summary
from server.models import Variants, Publications, VariantsAcmgRules, Classification, Reviews, AcmgRules, \
    ScientificMembers
from server.responses.internal_response import InternalResponse


def load_review_page_content(vus_id: str) -> Tuple[Dict, List[Dict], List[Dict], List[str]]:
    vus: Variants = db.session.query(Variants).get(vus_id)

    # TODO: remove extra info not used in front
    vus_publications: List[Publications] = vus.publications
    publications = [{"id": p.id, "title": p.title, "doi": p.doi} for p in vus_publications]

    vus_acmg_rules: List[VariantsAcmgRules] = vus.variants_acmg_rules
    acmg_rules = [{"id": r.acmg_rule_id, "name": r.rule_name.value} for r in vus_acmg_rules]

    variant_summary = get_variant_summary(vus)
    variant_summary['classification'] = vus.classification.value

    classifications = [Classification.BENIGN.value, Classification.LIKELY_BENIGN.value,
                       Classification.VUS.value, Classification.LIKELY_PATHOGENIC.value,
                       Classification.PATHOGENIC.value]

    return variant_summary, publications, acmg_rules, classifications


def save_review(vus_id: str, new_classification: str, reason: str, publication_ids: List[int],
                acmg_rule_ids: List[int]):
    review = Reviews(variant_id=vus_id, scientific_member_id=current_user.id, date_added=datetime.now(),
                     classification=new_classification.replace(" ", "_"), classification_reason=reason)

    vus: Variants = db.session.query(Variants).get(vus_id)
    publications: List[Publications] = db.session.query(Publications).filter(Publications.id.in_(publication_ids)).all()
    acmg_rules: List[AcmgRules] = db.session.query(AcmgRules).filter(AcmgRules.id.in_(acmg_rule_ids)).all()

    review.variant = vus
    review.publications = publications
    review.acmg_rules = acmg_rules

    vus.classification = new_classification.replace(" ", "_")

    db.session.add(review)

    try:
        # Commit the session to persist changes to the database
        db.session.commit()
        return InternalResponse({'isSuccess': True}, 200)
    except SQLAlchemyError as e:
        # Changes were rolled back due to an error
        db.session.rollback()

        current_app.logger.error(
            f'Rollback carried out since insertion of Classification Review entry for variant with id {vus_id} by user with id {current_user.id} in DB failed due to error: {e}')
        return InternalResponse({'isSuccess': False}, 500)


def get_all_reviews(vus_id: str):
    vus: Variants = db.session.query(Variants).get(vus_id)

    reviews: List[Reviews] = sorted(vus.reviews, key=lambda x: x.date_added, reverse=True)

    reviews_list = []
    for review in reviews:
        scientific_member: ScientificMembers = db.session.query(ScientificMembers).get(review.scientific_member_id)

        acmg_rules_list = []
        if review.acmg_rules:
            acmg_rules_list = [r.rule_name.value for r in review.acmg_rules]

        publications_list = []
        if review.publications:
            publications_list = [{"title": p.title, "doi": p.doi, "link": p.link} for p in review.publications]

        r = {"classification": review.classification.value, "reason": review.classification_reason,
             "dateAdded": review.date_added.strftime('%Y/%m/%d'), "scientificMemberName": scientific_member.name + " " +
                                                                     scientific_member.surname,
             "scientificMemberEmail": scientific_member.email, "acmgRules": acmg_rules_list,
             "publications": publications_list}

        reviews_list.append(r)

    variant_summary = get_variant_summary(vus)

    return reviews_list, variant_summary, vus.date_added.strftime('%Y/%m/%d')
