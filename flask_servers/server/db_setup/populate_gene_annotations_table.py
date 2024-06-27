import gzip
from typing import List, Tuple

from flask import current_app
from sqlalchemy import text
from sqlalchemy.exc import SQLAlchemyError

from server import db
from server.models import GeneAnnotations, GeneAttributes


def parse_attributes(gene_id: int, attributes: str) -> List[Tuple[int, str, str]]:
    separated_attributes = attributes.split('; ')

    attribute_list = []

    for attribute in separated_attributes:
        if attribute.strip() != '':
            attribute_split = attribute.split()
            attribute_list.append((gene_id, attribute_split[0], attribute_split[1].replace("\"", "").replace(';', '')))

    return attribute_list


def store_gtf_file_in_db():
    num_of_rows_added = 0

    with gzip.open('server/db_setup/Homo_sapiens.GRCh37.87.gtf.gz', 'rt') as gtf_file:
        for annotation_line in gtf_file:
            # skip comments
            if not annotation_line.startswith('#'):
                split_annotation_line = annotation_line.strip().split('\t')

                # if feature is a gene
                if split_annotation_line[2].lower() == 'gene':
                    num_of_rows_added += 1

                    # Replace '.' with None
                    updated_split_annotation_line = [None if info == '.' else info for info in split_annotation_line]

                    # Replace strand with 'POSITIVE' or 'NEGATIVE'
                    if updated_split_annotation_line[6] == '+':
                        updated_split_annotation_line[6] = 'POSITIVE'
                    else:
                        updated_split_annotation_line[6] = 'NEGATIVE'

                    if updated_split_annotation_line[5] is not None:
                        score = updated_split_annotation_line[5]
                    else:
                        score = "0"

                    try:
                        # Insert gene annotation
                        gene_annotation = GeneAnnotations(
                            seq_name=updated_split_annotation_line[0],
                            source=updated_split_annotation_line[1],
                            feature=updated_split_annotation_line[2],
                            start_location=int(updated_split_annotation_line[3]),
                            end_location=int(updated_split_annotation_line[4]),
                            score=float(score),
                            strand=updated_split_annotation_line[6],
                            frame=updated_split_annotation_line[7]
                        )
                        db.session.add(gene_annotation)
                        db.session.flush()  # Ensure the object is written to the database and obtain its id

                        gene_id = gene_annotation.id

                        # Parse and insert gene attributes
                        parsed_attributes = parse_attributes(gene_id, updated_split_annotation_line[8])
                        for attribute in parsed_attributes:
                            gene_attribute = GeneAttributes(
                                gene_id=attribute[0],
                                attribute_name=attribute[1],
                                attribute_value=attribute[2]
                            )
                            db.session.add(gene_attribute)

                        db.session.commit()

                    except SQLAlchemyError as e:
                        current_app.logger.error(f'Error {e}: on line {updated_split_annotation_line}')
                        db.session.rollback()


