from typing import List

import pandas as pd

from server import db
from server.models import Samples, SampleFiles, VariantsSamples, t_samples_phenotypes, Phenotypes


def retrieve_all_samples_from_db():
    samples_arr = []

    # retrieve all samples
    all_samples = db.session.query(Samples).all()

    # for each sample retrieve the file information and the variants it has
    for sample in all_samples:
        sample_file: SampleFiles = (
            db.session.query(SampleFiles).filter(SampleFiles.sample_file_id == sample.sample_file_id)
            .one_or_none())

        # Query with filter condition on sample_id
        phenotype_ontology_term_ids_res: List[str] = (db.session.query(t_samples_phenotypes.c.ontology_term_id)
                                                  .filter(t_samples_phenotypes.c.sample_id == sample.sample_id).all())

        phenotype_ontology_term_ids = [x[0] for x in phenotype_ontology_term_ids_res]

        sample_phenotypes: List[Phenotypes] = (db.session.query(Phenotypes)
                                        .filter(Phenotypes.ontology_term_id.in_(phenotype_ontology_term_ids)).all())

        phenotypes = []

        for p in sample_phenotypes:
            sample_phenotype = {'ontologyId': p.ontology_term_id, 'name': p.term_name}
            phenotypes.append(sample_phenotype)

        variants_samples: List[VariantsSamples] = db.session.query(VariantsSamples).filter(
            VariantsSamples.sample_id == sample.sample_id).all()

        variants = []

        for v in variants_samples:
            variant_sample = {'variantId': v.variant_id, 'genotype': v.genotype.value}
            variants.append(variant_sample)

        samples_arr.append({'sampleId': sample.sample_id, 'phenotype': phenotypes,
                            'genomeVersion': sample.genome_version, 'fileUploadName': sample_file.filename,
                            'dateOfFileUpload': str(sample_file.date_uploaded.date()), 'variants': variants})

    return samples_arr
