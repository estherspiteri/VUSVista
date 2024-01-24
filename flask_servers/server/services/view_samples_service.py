from typing import List

import pandas as pd

from server import db
from server.models import Samples, SampleFiles, VariantsSamples


def retrieve_all_samples_from_db():
    samples_arr = []

    # retrieve all samples
    all_samples = db.session.query(Samples).all()

    # for each sample retrieve the file information and the variants it has
    for sample in all_samples:
        sample_file: SampleFiles = (db.session.query(SampleFiles).filter(SampleFiles.sample_file_id == sample.sample_file_id)
                       .one_or_none())

        variants_samples: List[VariantsSamples] = db.session.query(VariantsSamples).filter(VariantsSamples.sample_id == sample.sample_id).all()

        variants = []

        for v in variants_samples:
            variant_sample = {'variantId': v.variant_id, 'genotype': v.genotype.value}
            variants.append(variant_sample)

        samples_arr.append({'sampleId': sample.sample_id, 'phenotype': sample.phenotype,
                            'genomeVersion': sample.genome_version, 'fileUploadName': sample_file.filename,
                            'dateOfFileUpload': str(sample_file.date_uploaded.date()), 'variants': variants})

    return samples_arr

