import math
from datetime import datetime

from flask_login import current_user

from server import db
from server.models import VariantsSamplesUploads, FileUploads, ManualUploads, VariantsSamples, VariantHgvs


def add_variant_sample_to_db(variant_id: int, sample_id: str, hgvs: str, genotype: str, consequence: str):
    # check if variants_samples entry already exists
    existing_variants_samples: VariantsSamples = db.session.query(VariantsSamples).filter(
        VariantsSamples.variant_id == variant_id, VariantsSamples.sample_id == sample_id).one_or_none()

    # if variants_samples entry does not exist, add new entry
    if existing_variants_samples is None:
        hgvs_id = None
        if not (isinstance(hgvs, float) and math.isnan(hgvs)):
            # check if HGVS exists
            variant_hgvs: VariantHgvs = db.session.query(VariantHgvs).filter(VariantHgvs.variant_id == variant_id,
                                                                             VariantHgvs.hgvs == hgvs).one_or_none()

            if variant_hgvs is None and hgvs is not None:
                # create new HGVS entry
                variant_hgvs = VariantHgvs(variant_id=variant_id, hgvs=hgvs, is_updated=False)
                db.session.add(variant_hgvs)
                db.session.flush()
                hgvs_id = variant_hgvs.id

        new_variants_samples = VariantsSamples(variant_id=variant_id, sample_id=sample_id,
                                               variant_hgvs_id=hgvs_id, genotype=genotype.upper(), consequence=consequence)
        db.session.add(new_variants_samples)


def store_upload_details_for_variant_sample(file_upload: FileUploads | None, is_file_upload: bool, sample_id: str,
                                            variant_id: int, user_id: int):
    # create sample upload entry
    if is_file_upload:
        upload_type = 'file'
    else:
        upload_type = 'manual'

    new_variants_sample_upload = VariantsSamplesUploads(date_uploaded=datetime.now(),
                                                        scientific_member_id=user_id,
                                                        upload_type=upload_type, sample_id=sample_id,
                                                        variant_id=variant_id)

    if is_file_upload and file_upload is not None:
        new_variants_sample_upload.file_upload = file_upload

    db.session.add(new_variants_sample_upload)
    db.session.flush()

    if not is_file_upload:
        # create manual uploads entry
        new_manual_upload = ManualUploads(variants_samples_uploads_manual_id=new_variants_sample_upload.id)
        db.session.add(new_manual_upload)