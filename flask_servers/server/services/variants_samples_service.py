from datetime import datetime

from flask_login import current_user

from server import db
from server.models import VariantsSamplesUploads, FileUploads, ManualUploads


def store_upload_details_for_variant_sample(file_upload: FileUploads | None, is_file_upload: bool, sample_id: str,
                                            variant_id: int):
    # create sample upload entry
    if is_file_upload:
        upload_type = 'file'
    else:
        upload_type = 'manual'

    new_variants_sample_upload = VariantsSamplesUploads(date_uploaded=datetime.now(),
                                                        scientific_member_id=current_user.id,
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