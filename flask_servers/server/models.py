## https://pypi.org/project/sqlacodegen-v2/
from typing import List

from enum import Enum

from sqlalchemy import Boolean, CHAR, Column, Date, DateTime, Double, Enum as EnumSQL, ForeignKeyConstraint, Identity, \
    Integer, PrimaryKeyConstraint, String, Table, Text, UniqueConstraint
from sqlalchemy.orm import declarative_base, mapped_column, relationship
from sqlalchemy.orm.base import Mapped

Base = declarative_base()
metadata = Base.metadata


# class ExternalRefDbType(Enum):
#     DBSNP = 'DBSNP'
#     CLINVAR = 'CLINVAR'


class VariantType(Enum):
    SNV = 'SNV'
    MNV = 'MNV'
    INDEL = 'INDEL'


class Consequence(Enum):
    MISSENSE = 'MISSENSE'
    NONSENSE = 'NONSENSE'
    INSERTION = 'INSERTION'
    DELETION = 'DELETION'
    FRAMESHIFT = 'FRAMESHIFT'
    DUPLICATION = 'DUPLICATION'


class Strand(Enum):
    POSITIVE = 'POSITIVE'
    NEGATIVE = 'NEGATIVE'


class Genotype(Enum):
    HOMOZYGOUS = 'HOMOZYGOUS'
    HETEROZYGOUS = 'HETEROZYGOUS'


class ACMGRule(Enum):
    PS2 = 'PS2'
    PM3 = 'PM3'
    PM6 = 'PM6'
    PP1 = 'PP1'
    PP4 = 'PP4'
    BS4 = 'BS4'
    BP2 = 'BP2'
    PS3 = 'PS3'
    PP5 = 'PP5'
    BP6 = 'BP6'


class ACMGStrength(Enum):
    SUPPORTING = 'SUPPORTING'
    MODERATE = 'MODERATE'
    STRONG = 'STRONG'
    VERY_STRONG = 'VERY_STRONG'


class Classification(Enum):
    PATHOGENIC = 'PATHOGENIC'
    LIKELY_PATHOGENIC = 'LIKELY_PATHOGENIC'
    VUS = 'VUS'
    UNCERTAIN_SIGNIFICANCE = 'UNCERTAIN_SIGNIFICANCE'
    UNCLASSIFIED = 'UNCLASSIFIED'
    LIKELY_VUS = 'LIKELY_VUS'
    LIKELY_BENIGN = 'LIKELY_BENIGN'
    BENIGN = 'BENIGN'


class ReviewStatus(Enum):
    IN_PROGRESS = 'IN_PROGRESS'
    COMPLETE = 'COMPLETE'


class AcmgRules(Base):
    __tablename__ = 'acmg_rules'
    __table_args__ = (
        PrimaryKeyConstraint('rule_name', name='acmg_rules_pkey'),
    )

    rule_name = mapped_column(EnumSQL(ACMGRule, name='acmg_rule'))
    description = mapped_column(Text, nullable=False)
    default_strength = mapped_column(EnumSQL(ACMGStrength, name='acmg_strength'), nullable=False)
    requires_lab_verification = mapped_column(Boolean, nullable=False)

    review: Mapped['Reviews'] = relationship('Reviews', secondary='reviews_acmg_rules', back_populates='acmg_rules')
    samples_variants_acmg_rules: Mapped[List['SamplesVariantsAcmgRules']] = relationship('SamplesVariantsAcmgRules', uselist=True, back_populates='acmg_rules')


class ExternalReferences(Base):
    __tablename__ = 'external_references'
    __table_args__ = (
        ForeignKeyConstraint(['variant_id'], ['variants.variant_id'], name='fk_variants'),
        PrimaryKeyConstraint('external_references_id', name='external_references_pkey')
    )

    external_references_id = mapped_column(Integer, Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False, cache=1))
    db_type = mapped_column(Text, nullable=False)
    variant_id = mapped_column(Integer, nullable=False)
    error_msg = mapped_column(Text)

    variant: Mapped['Variants'] = relationship('Variants', back_populates='external_references')
    clinvar: Mapped['Clinvar'] = relationship('Clinvar', uselist=False, back_populates='external_clinvar')
    db_snp: Mapped['DbSnp'] = relationship('DbSnp', uselist=False, back_populates='external_db_snp')


class Clinvar(Base):
    __tablename__ = 'clinvar'
    __table_args__ = (
        ForeignKeyConstraint(['external_clinvar_id'], ['external_references.external_references_id'], name='fk_external_references'),
        PrimaryKeyConstraint('clinvar_id', name='clinvar_pkey'),
        UniqueConstraint('external_clinvar_id', name='clinvar_external_clinvar_id_key')
    )

    clinvar_id = mapped_column(Text)
    external_clinvar_id = mapped_column(Integer, nullable=False)
    canonical_spdi = mapped_column(Text, nullable=False)
    classification = mapped_column(Text)
    last_evaluated = mapped_column(DateTime)
    review_status = mapped_column(Text)

    external_clinvar: Mapped['ExternalReferences'] = relationship('ExternalReferences', back_populates='clinvar')


class DbSnp(Base):
    __tablename__ = 'db_snp'
    __table_args__ = (
        ForeignKeyConstraint(['external_db_snp_id'], ['external_references.external_references_id'], name='fk_external_references'),
        PrimaryKeyConstraint('db_snp_id', name='db_snp_pkey'),
        UniqueConstraint('external_db_snp_id', name='db_snp_external_db_snp_id_key')
    )

    db_snp_id = mapped_column(String(15))
    external_db_snp_id = mapped_column(Integer, nullable=False)

    external_db_snp: Mapped['ExternalReferences'] = relationship('ExternalReferences', back_populates='db_snp')


class GeneAnnotations(Base):
    __tablename__ = 'gene_annotations'
    __table_args__ = (
        PrimaryKeyConstraint('gene_id', name='gene_annotations_pkey'),
    )

    gene_id = mapped_column(Integer, Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False, cache=1))
    seq_name = mapped_column(Text, nullable=False)
    source = mapped_column(Text, nullable=False)
    feature = mapped_column(Text, nullable=False)
    start_location = mapped_column(Integer, nullable=False)
    end_location = mapped_column(Integer, nullable=False)
    score = mapped_column(Double(53))
    strand = mapped_column(EnumSQL(Strand, name='strand'))
    frame = mapped_column(CHAR(1))

    gene_attributes: Mapped[List['GeneAttributes']] = relationship('GeneAttributes', uselist=True, back_populates='gene')
    variants: Mapped[List['Variants']] = relationship('Variants', uselist=True, back_populates='gene')


class Publications(Base):
    __tablename__ = 'publications'
    __table_args__ = (
        PrimaryKeyConstraint('pmid', name='publications_pkey'),
    )

    pmid = mapped_column(Integer)
    title = mapped_column(Text, nullable=False)
    match_in_sup_material = mapped_column(Boolean, nullable=False)
    date_published = mapped_column(Date, nullable=False)
    abstract = mapped_column(Text)

    variant: Mapped['Variants'] = relationship('Variants', secondary='variants_publications', back_populates='publications')
    review: Mapped['Reviews'] = relationship('Reviews', secondary='reviews_publications', back_populates='publications')


class SampleFiles(Base):
    __tablename__ = 'sample_files'
    __table_args__ = (
        PrimaryKeyConstraint('sample_file_id', name='sample_files_pkey'),
    )

    sample_file_id = mapped_column(Integer, Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False, cache=1))
    date_uploaded = mapped_column(DateTime, nullable=False)
    filename = mapped_column(Text, nullable=False)

    samples: Mapped[List['Samples']] = relationship('Samples', uselist=True, back_populates='sample_file')


class ScientificMembers(Base):
    __tablename__ = 'scientific_members'
    __table_args__ = (
        PrimaryKeyConstraint('scientific_member_id', name='scientific_members_pkey'),
    )

    scientific_member_id = mapped_column(Integer, Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False, cache=1))
    name = mapped_column(Text, nullable=False)
    surname = mapped_column(Text, nullable=False)

    reviews: Mapped[List['Reviews']] = relationship('Reviews', uselist=True, back_populates='scientific_member')


class GeneAttributes(Base):
    __tablename__ = 'gene_attributes'
    __table_args__ = (
        ForeignKeyConstraint(['gene_id'], ['gene_annotations.gene_id'], name='fk_gene_annotations'),
        PrimaryKeyConstraint('gene_id', 'attribute_name', name='gene_attributes_pkey')
    )

    gene_id = mapped_column(Integer, nullable=False)
    attribute_name = mapped_column(Text, nullable=False)
    attribute_value = mapped_column(Text, nullable=False)

    gene: Mapped['GeneAnnotations'] = relationship('GeneAnnotations', back_populates='gene_attributes')


class Samples(Base):
    __tablename__ = 'samples'
    __table_args__ = (
        ForeignKeyConstraint(['sample_file_id'], ['sample_files.sample_file_id'], name='fk_sample_files'),
        PrimaryKeyConstraint('sample_id', name='samples_pkey')
    )

    sample_id = mapped_column(Integer, Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False, cache=1))
    sample_file_id = mapped_column(Integer, nullable=False)
    phenotype = mapped_column(Text)
    genome_version = mapped_column(String(20))

    sample_file: Mapped['SampleFiles'] = relationship('SampleFiles', back_populates='samples')
    samples_variants_acmg_rules: Mapped[List['SamplesVariantsAcmgRules']] = relationship('SamplesVariantsAcmgRules', uselist=True, back_populates='sample')
    variants_samples: Mapped[List['VariantsSamples']] = relationship('VariantsSamples', uselist=True, back_populates='sample')


class Variants(Base):
    __tablename__ = 'variants'
    __table_args__ = (
        ForeignKeyConstraint(['gene_id'], ['gene_annotations.gene_id'], name='fk_gene_annotations'),
        PrimaryKeyConstraint('variant_id', name='variants_pkey')
    )

    variant_id = mapped_column(Integer, Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False, cache=1))
    chromosome = mapped_column(String(2), nullable=False)
    chromosome_position = mapped_column(Text, nullable=False)
    variant_type = mapped_column(EnumSQL(VariantType, name='variant_type'), nullable=False)
    ref = mapped_column(Text, nullable=False)
    classification = mapped_column(EnumSQL(Classification, name='classification'), nullable=False)
    gene_id = mapped_column(Integer, nullable=False)
    gene_name = mapped_column(Text, nullable=False)
    alt = mapped_column(Text)
    consequences = mapped_column(EnumSQL(Consequence, name='consequence'))

    publications: Mapped['Publications'] = relationship('Publications', secondary='variants_publications', back_populates='variant')
    gene: Mapped['GeneAnnotations'] = relationship('GeneAnnotations', back_populates='variants')
    external_references: Mapped[List['ExternalReferences']] = relationship('ExternalReferences', uselist=True, back_populates='variant')
    reviews: Mapped[List['Reviews']] = relationship('Reviews', uselist=True, back_populates='variant')
    samples_variants_acmg_rules: Mapped[List['SamplesVariantsAcmgRules']] = relationship('SamplesVariantsAcmgRules', uselist=True, back_populates='variant')
    variants_samples: Mapped[List['VariantsSamples']] = relationship('VariantsSamples', uselist=True, back_populates='variant')


class Reviews(Base):
    __tablename__ = 'reviews'
    __table_args__ = (
        ForeignKeyConstraint(['scientific_member_id'], ['scientific_members.scientific_member_id'], name='fk_scientific_members'),
        ForeignKeyConstraint(['variant_id'], ['variants.variant_id'], name='fk_variants'),
        PrimaryKeyConstraint('review_id', name='reviews_pkey')
    )

    review_id = mapped_column(Integer, Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False, cache=1))
    variant_id = mapped_column(Integer, nullable=False)
    scientific_member_id = mapped_column(Integer, nullable=False)
    date_added = mapped_column(DateTime, nullable=False)
    classification = mapped_column(EnumSQL(Classification, name='classification'), nullable=False)
    review_status = mapped_column(EnumSQL(ReviewStatus, name='review_status'))
    classification_reason = mapped_column(Text)


    publications: Mapped['Publications'] = relationship('Publications', secondary='reviews_publications', back_populates='review')
    scientific_member: Mapped['ScientificMembers'] = relationship('ScientificMembers', back_populates='reviews')
    variant: Mapped['Variants'] = relationship('Variants', back_populates='reviews')
    acmg_rules: Mapped['AcmgRules'] = relationship('AcmgRules', secondary='reviews_acmg_rules', back_populates='review')


class SamplesVariantsAcmgRules(Base):
    __tablename__ = 'samples_variants_acmg_rules'
    __table_args__ = (
        ForeignKeyConstraint(['rule_name'], ['acmg_rules.rule_name'], name='fk_acmg_rules'),
        ForeignKeyConstraint(['sample_id'], ['samples.sample_id'], name='fk_samples'),
        ForeignKeyConstraint(['variant_id'], ['variants.variant_id'], name='fk_variants'),
        PrimaryKeyConstraint('variant_id', 'sample_id', 'rule_name', name='samples_variants_acmg_rules_pkey')
    )

    variant_id = mapped_column(Integer, nullable=False)
    sample_id = mapped_column(Integer, nullable=False)
    rule_name = mapped_column(EnumSQL(ACMGRule, name='acmg_rule'), nullable=False)

    acmg_rules: Mapped['AcmgRules'] = relationship('AcmgRules', back_populates='samples_variants_acmg_rules')
    sample: Mapped['Samples'] = relationship('Samples', back_populates='samples_variants_acmg_rules')
    variant: Mapped['Variants'] = relationship('Variants', back_populates='samples_variants_acmg_rules')


t_variants_publications = Table(
    'variants_publications', metadata,
    Column('variant_id', Integer, nullable=False),
    Column('pmid', Integer, nullable=False),
    ForeignKeyConstraint(['pmid'], ['publications.pmid'], name='fk_publications'),
    ForeignKeyConstraint(['variant_id'], ['variants.variant_id'], name='fk_variants'),
    PrimaryKeyConstraint('variant_id', 'pmid', name='variants_publications_pkey')
)


class VariantsSamples(Base):
    __tablename__ = 'variants_samples'
    __table_args__ = (
        ForeignKeyConstraint(['sample_id'], ['samples.sample_id'], name='fk_samples'),
        ForeignKeyConstraint(['variant_id'], ['variants.variant_id'], name='fk_variants'),
        PrimaryKeyConstraint('variant_id', 'sample_id', name='variants_samples_pkey')
    )

    variant_id = mapped_column(Integer, nullable=False)
    sample_id = mapped_column(Integer, nullable=False)
    genotype = mapped_column(EnumSQL(Genotype, name='genotype'), nullable=False)

    sample: Mapped['Samples'] = relationship('Samples', back_populates='variants_samples')
    variant: Mapped['Variants'] = relationship('Variants', back_populates='variants_samples')

t_reviews_acmg_rules = Table(
    'reviews_acmg_rules', metadata,
    Column('review_id', Integer, nullable=False),
    Column('rule_name', EnumSQL(ACMGRule, name='acmg_rule'), nullable=False),
    ForeignKeyConstraint(['review_id'], ['reviews.review_id'], name='fk_reviews'),
    ForeignKeyConstraint(['rule_name'], ['acmg_rules.rule_name'], name='fk_acmg_rules'),
    PrimaryKeyConstraint('review_id', 'rule_name', name='reviews_acmg_rules_pkey')
)


t_reviews_publications = Table(
    'reviews_publications', metadata,
    Column('review_id', Integer, nullable=False),
    Column('pmid', Integer, nullable=False),
    ForeignKeyConstraint(['pmid'], ['publications.pmid'], name='fk_publications'),
    ForeignKeyConstraint(['review_id'], ['reviews.review_id'], name='fk_reviews'),
    PrimaryKeyConstraint('review_id', 'pmid', name='reviews_publications_pkey')
)

