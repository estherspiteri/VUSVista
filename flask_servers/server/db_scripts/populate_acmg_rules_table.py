import psycopg2

from server.models import ACMGStrength


def store_acmg_rules_in_db(cursor) -> bool:
   # rule_name, description, default_strength, requires_lab_verification
    acmg_rules = [
        ('PS2',
        'De novo (both maternity and paternity confirmed) in a patient with the disease and no family history.',
        str(ACMGStrength.STRONG.value),
        True
         ),
        ('PM3',
        'For recessive disorders, detected in trans with a pathogenic variant.',
        str(ACMGStrength.MODERATE.value),
        True
        ),
        ('PM6',
        'Assumed de novo, but without confirmation of paternity and maternity.',
        str(ACMGStrength.MODERATE.value),
        True
         ),
       ('PP1',
        'Cosegregation with disease in multiple affected family members in a gene definitively known to cause the '
        'disease.',
        str(ACMGStrength.SUPPORTING.value),
        True
        ),
        ('PP4',
        'Patient’s phenotype or family history is highly specific for a disease with a single genetic etiology.',
        str(ACMGStrength.SUPPORTING.value),
        True
         ),
        ('BS4',
        'Lack of segregation in affected members of a family.',
        str(ACMGStrength.STRONG.value),
        True
         ),
        ('BP2',
        'Observed in trans with a pathogenic variant for a fully penetrant dominant gene/disorder or observed in cis '
        'with a pathogenic variant in any inheritance pattern.',
        str(ACMGStrength.SUPPORTING.value),
        True
         ),
       # associated to literature
        ('PS3',
        'Well-established in vitro or in vivo functional studies supportive of a damaging effect on the gene or gene '
        'product.',
        str(ACMGStrength.STRONG.value),
        False
         ),
        ('PP5',
        'Reputable source recently reports variant as pathogenic, but the evidence is not available to the laboratory '
        'to perform an independent evaluation.',
        str(ACMGStrength.SUPPORTING.value),
        False
         ),
        ('BP6',
        'Reputable source recently reports variant as benign, but the evidence is not available to the laboratory to '
        'perform an independent evaluation.',
        str(ACMGStrength.SUPPORTING.value),
        False
         )]

    for rule in acmg_rules:
        try:
            sql_string = ("INSERT INTO acmg_rules(rule_name, description, default_strength, requires_lab_verification) "
                          "VALUES(%s,%s,%s,%s);")

            # inserting gene annotations in db
            cursor.execute(sql_string,  rule)

        except BaseException as e:
            print(f'Error {e}: whilst inputting ACMG rule {rule[0]} into database')
            return False

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

    # emptying acmg_rules table
    cursor.execute("TRUNCATE acmg_rules CASCADE")

    is_success = store_acmg_rules_in_db(cursor)

    if not is_success:
        # rollback changes
        conn.rollback()
    else:
        # committing changes
        conn.commit()

    # closing connection
    conn.close()


main()