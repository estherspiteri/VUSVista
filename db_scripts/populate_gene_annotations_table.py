import gzip
from typing import List, Tuple

import psycopg2


def parse_attributes(gene_id: int, attributes: str) -> List[Tuple[str, str, str]]:
    separated_attributes = attributes.split('; ')

    attribute_list = []

    for attribute in separated_attributes:
        attribute_split = attribute.split()
        attribute_list.append((str(gene_id), attribute_split[0],
                               attribute_split[1].replace("\"", "").replace(';', '')))

    return attribute_list


def store_gtf_file_in_db(zip_filename: str, cursor) -> bool:
    num_of_rows_added = 0

    with gzip.open(zip_filename, 'rt') as gtf_file:
        for annotation_line in gtf_file:
            # skip comments
            if not annotation_line.startswith('#'):
                split_annotation_line = annotation_line.strip().split('\t')

                # if feature is a gene
                if split_annotation_line[2].lower() == 'gene':
                    num_of_rows_added += 1
                    updated_split_annotation_line = []

                    # where value is not set, replace '.' with None
                    for info in split_annotation_line:
                        if info == '.':
                            updated_split_annotation_line.append(None)
                        else:
                            updated_split_annotation_line.append(info)

                    # replace strand with 'positive' or 'negative'
                    if updated_split_annotation_line[6] == '+':
                        updated_split_annotation_line[6] = 'POSITIVE'
                    else:
                        updated_split_annotation_line[6] = 'NEGATIVE'

                    try:
                        sql_string = ("INSERT INTO gene_annotations(seq_name, source, feature, start_location, "
                                      "end_location, score, strand, frame) VALUES(%s,%s,%s,%s,%s,%s,%s,%s) "
                                      "RETURNING id;")

                        # inserting gene annotations in db
                        cursor.execute(sql_string, updated_split_annotation_line[0:8])

                        gene_id = cursor.fetchone()[0]

                        parsed_attributes = parse_attributes(gene_id, updated_split_annotation_line[8])

                        # inserting gene attributes in db
                        args_str = ','.join(cursor.mogrify("(%s,%s,%s)", a).decode('utf-8') for a in parsed_attributes)
                        cursor.execute("INSERT INTO gene_attributes VALUES " + args_str)

                    except BaseException as e:
                        print(f'Error {e}: on line {updated_split_annotation_line}')
                        return False

    print(f'Number of rows added: {num_of_rows_added}')
    return True


def main():
    conn = psycopg2.connect(
        host="localhost",
        database="vus-app-db",
        user="postgres",
        password="21641" #TODO: store securely
        )

    # creating a cursor
    cursor = conn.cursor()

    # emptying gene_annotation table
    cursor.execute("TRUNCATE gene_annotations CASCADE")

    is_success = store_gtf_file_in_db('Homo_sapiens.GRCh37.87.gtf.gz', cursor)

    if not is_success:
        # rollback changes
        conn.rollback()
    else:
        # committing changes
        conn.commit()

    # closing connection
    conn.close()


main()