from typing import List, Dict

import pandas as pd

from server import db
from server.helpers.data_helper import get_variant_summary
from server.models import Samples, VariantsSamples, t_samples_phenotypes, Phenotypes, \
    Variants, VariantsSamplesUploads


def get_sample_info_from_db(sample: Samples) -> Dict:
    # Query with filter condition on sample_id
    phenotype_ontology_term_ids_res: List[str] = (db.session.query(t_samples_phenotypes.c.ontology_term_id)
                                                  .filter(t_samples_phenotypes.c.sample_id == sample.id).all())

    phenotype_ontology_term_ids = [x[0] for x in phenotype_ontology_term_ids_res]

    sample_phenotypes: List[Phenotypes] = (db.session.query(Phenotypes)
                                           .filter(Phenotypes.ontology_term_id.in_(phenotype_ontology_term_ids)).all())

    phenotypes = []

    for p in sample_phenotypes:
        sample_phenotype = {'ontologyId': p.ontology_term_id, 'name': p.term_name}
        phenotypes.append(sample_phenotype)

    variants_samples: List[VariantsSamples] = db.session.query(VariantsSamples).filter(
        VariantsSamples.sample_id == sample.id).all()

    variants = []

    for v_s in variants_samples:
        variant_details: Variants = db.session.query(Variants).filter(Variants.id == v_s.variant_id).first()

        variant_summary = get_variant_summary(variant_details)

        # get upload files
        files: List[VariantsSamplesUploads] = [u for u in v_s.variants_samples_uploads if u.upload_type == 'file']

        file_dicts = [{'filename': u.file_upload.filename, 'dateOfFileUpload': str(u.date_uploaded.date())} for u in
                      files]

        variant_sample = {'variantId': v_s.variant_id, 'variant': variant_summary, 'genotype': v_s.genotype.value,
                          'files': file_dicts}
        variants.append(variant_sample)

    return {'sampleId': sample.id, 'phenotype': phenotypes, 'genomeVersion': sample.genome_version,
            'variants': variants}


def retrieve_all_samples_from_db():
    samples_arr = []

    # retrieve all samples
    all_samples = db.session.query(Samples).all()

    # for each sample retrieve the number of variants it has
    for sample in all_samples:
        variants_count: List[VariantsSamples] = db.session.query(VariantsSamples).filter(
            VariantsSamples.sample_id == sample.id).count()

        samples_arr.append({'sampleId': sample.id, 'numOfVariants': variants_count})

    return samples_arr


def retrieve_sample_from_db(sample_id: str):
    # retrieve sample
    sample = db.session.query(Samples).filter(Samples.id == sample_id).first()
    sample_dict = get_sample_info_from_db(sample)

    return sample_dict
