import psycopg2

from server.models import ACMGStrength


def store_acmg_rules_in_db(cursor) -> bool:
   # rule_name, description, default_strength, requires_lab_verification
    acmg_rules = [
        # associated to family history/disease knowledge
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
         ),
        # remaining rules
        ('PVS1',
         'Null variant (nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single or multiexon '
         'deletion) in a gene where LOF is a known mechanism of disease.',
         str(ACMGStrength.VERY_STRONG.value),
         False
         ),
        ('PS1',
         'Same amino acid change as a previously established pathogenic variant regardless of nucleotide change.',
         str(ACMGStrength.STRONG.value),
         False
         ),
        ('BS3',
         'Well-established in vitro or in vivo functional studies show no damaging effect on protein function or '
         'splicing. ',
         str(ACMGStrength.STRONG.value),
         False
         ),
        ('PM1',
         'Located in a mutational hot spot and/or critical and well-established functional domain (e.g., active site '
         'of an enzyme) without benign variation.',
         str(ACMGStrength.MODERATE.value),
         False
         ),
        ('BP3',
         'In-frame deletions/insertions in a repetitive region without a known function.',
         str(ACMGStrength.SUPPORTING.value),
         False
         ),
        ('PM2',
         'Absent from controls (or at extremely low frequency if recessive) in Exome Sequencing Project, 1000 Genomes '
         'Project, or Exome Aggregation Consortium.',
         str(ACMGStrength.MODERATE.value),
         False
         ),
        ('PM4',
         'Protein length changes as a result of in-frame deletions/insertions in a non-repeat region or stop-loss '
         'variants.',
         str(ACMGStrength.MODERATE.value),
         False
         ),
        ('PM5',
         'Novel missense change at an amino acid residue where a different missense change determined to be pathogenic '
         'has been seen before. ',
         str(ACMGStrength.MODERATE.value),
         False
         ),
        ('PP2',
         'Missense variant in a gene that has a low rate of benign missense variation and in which missense variants '
         'are a common mechanism of disease.',
         str(ACMGStrength.SUPPORTING.value),
         False
         ),
        ('BP1',
         'Missense variant in a gene for which primarily truncating variants are known to cause disease.',
         str(ACMGStrength.SUPPORTING.value),
         False
         ),
        ('PP3',
         'Multiple lines of computational evidence support a deleterious effect on the gene or gene product('
         'conservation, evolutionary, splicing impact, etc.)',
         str(ACMGStrength.SUPPORTING.value),
         False
         ),
        ('BP4',
         'Multiple lines of computational evidence suggest no impact on gene or gene product (conservation, '
         'evolutionary, splicing impact, etc.)',
         str(ACMGStrength.SUPPORTING.value),
         False
         ),
        ('BA1',
         'Allele frequency is >5% in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium.',
         str(ACMGStrength.STAND_ALONE.value),
         False
         ),
        ('BS1',
         'Allele frequency is greater than expected for disorder.',
         str(ACMGStrength.STRONG.value),
         False
         ),
        ('BS2',
         'Observed in a healthy adult individual for a recessive (homozygous), dominant (heterozygous), or X-linked '
         '(hemizygous) disorder, with full penetrance expected at an early age.',
         str(ACMGStrength.STRONG.value),
         False
         ),
        ('BP7',
         'A synonymous (silent) variant for which splicing prediction algorithms predict no impact to the splice '
         'consensus sequence nor the creation of a new splice site AND the nucleotide is not highly conserved. ',
         str(ACMGStrength.SUPPORTING.value),
         False
         ),
        ('PS4',
         'The prevalence of the variant in affected individuals is significantly increased compared with the '
         'prevalence in controls.',
         str(ACMGStrength.STRONG.value),
         False
         ),
        ('BP5',
         'Variant found in a case with an alternate molecular basis for disease.',
         str(ACMGStrength.SUPPORTING.value),
         True
         ),
    ]

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