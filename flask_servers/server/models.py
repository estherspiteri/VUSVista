## https://pypi.org/project/sqlacodegen-v2/
# cmd: sqlacodegen_v2 postgresql://postgres:21641@localhost:5432/vus-app-db
from typing import List

from enum import Enum

from sqlalchemy import Boolean, CHAR, Column, Date, DateTime, Double, Enum as EnumSQL, ForeignKeyConstraint, Identity, \
    Integer, PrimaryKeyConstraint, String, Table, Text, UniqueConstraint
from sqlalchemy.orm import declarative_base, mapped_column, relationship
from sqlalchemy.orm.base import Mapped, Optional
from flask_login import UserMixin

from server import db

Base = declarative_base()
Base.query = db.session.query_property()
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
    PVS1 = 'PVS1'
    PS1 = 'PS1'
    BS3 = 'BS3'
    PM1 = 'PM1'
    BP3 = 'BP3'
    PM2 = 'PM2'
    PM4 = 'PM4'
    PM5 = 'PM5'
    PP2 = 'PP2'
    BP1 = 'BP1'
    PP3 = 'PP3'
    BP4 = 'BP4'
    BA1 = 'BA1'
    BS1 = 'BS1'
    BS2 = 'BS2'
    BP7 = 'BP7'
    PS4 = 'PS4'
    BP5 = 'BP5'


class ACMGStrength(Enum):
    SUPPORTING = 'SUPPORTING'
    MODERATE = 'MODERATE'
    STRONG = 'STRONG'
    VERY_STRONG = 'VERY_STRONG'
    STAND_ALONE = 'STAND_ALONE'


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
        PrimaryKeyConstraint('id', name='acmg_rules_pkey'),
    )

    id = mapped_column(Integer, Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False, cache=1))
    rule_name = mapped_column(EnumSQL(ACMGRule, name='acmg_rule'))
    description = mapped_column(Text, nullable=False)
    default_strength = mapped_column(EnumSQL(ACMGStrength, name='acmg_strength'), nullable=False)
    requires_lab_verification = mapped_column(Boolean, nullable=False)

    review: Mapped['Reviews'] = relationship('Reviews', secondary='reviews_acmg_rules', back_populates='acmg_rules')
    variants_acmg_rules: Mapped[List['VariantsAcmgRules']] = relationship('VariantsAcmgRules', uselist=True,
                                                                          back_populates='acmg_rule')


class ExternalReferences(Base):
    __tablename__ = 'external_references'
    __table_args__ = (
        ForeignKeyConstraint(['variant_id'], ['variants.id'], name='fk_variants'),
        PrimaryKeyConstraint('id', name='external_references_pkey')
    )

    id = mapped_column(Integer,
                       Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False,
                                cache=1))
    db_type = mapped_column(Text, nullable=False)
    variant_id = mapped_column(Integer, nullable=False)
    error_msg = mapped_column(Text)

    variant: Mapped['Variants'] = relationship('Variants', back_populates='external_references')
    clinvar: Mapped['Clinvar'] = relationship('Clinvar', uselist=False, back_populates='external_clinvar')
    db_snp: Mapped['DbSnp'] = relationship('DbSnp', uselist=False, back_populates='external_db_snp')


class ClinvarUpdates(Base):
    __tablename__ = 'clinvar_updates'
    __table_args__ = (
        PrimaryKeyConstraint('id', name='clinvar_updates_pkey'),
    )

    id = mapped_column(Integer, Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False, cache=1))
    classification = mapped_column(Text)
    review_status = mapped_column(Text)
    last_evaluated = mapped_column(DateTime)

    clinvar_eval_dates: Mapped['ClinvarEvalDates'] = relationship('ClinvarEvalDates', back_populates='clinvar_update')


class ClinvarEvalDates(Base):
    __tablename__ = 'clinvar_eval_dates'
    __table_args__ = (
        ForeignKeyConstraint(['clinvar_id'], ['clinvar.id'], name='fk_clinvar'),
        ForeignKeyConstraint(['clinvar_update_id'], ['clinvar_updates.id'], name='fk_clinvar_updates'),
        PrimaryKeyConstraint('id', name='clinvar_eval_dates_pkey')
    )

    id = mapped_column(Integer, Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False, cache=1))
    clinvar_id = mapped_column(Integer, nullable=False)
    clinvar_update_id = mapped_column(Integer)
    eval_date = mapped_column(DateTime)

    clinvar: Mapped['Clinvar'] = relationship('Clinvar', back_populates='clinvar_eval_dates')
    clinvar_update: Mapped[Optional['ClinvarUpdates']] = relationship('ClinvarUpdates', back_populates='clinvar_eval_dates')


class Clinvar(Base):
    __tablename__ = 'clinvar'
    __table_args__ = (
        ForeignKeyConstraint(['external_clinvar_id'], ['external_references.id'], name='fk_external_references'),
        PrimaryKeyConstraint('id', name='clinvar_pkey'),
        UniqueConstraint('external_clinvar_id', name='clinvar_external_clinvar_id_key')
    )

    id = mapped_column(Integer, Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False, cache=1))
    uid = mapped_column(Text, nullable=False)
    external_clinvar_id = mapped_column(Integer, nullable=False)
    canonical_spdi = mapped_column(Text, nullable=False)

    external_clinvar: Mapped['ExternalReferences'] = relationship('ExternalReferences', back_populates='clinvar')
    clinvar_eval_dates: Mapped[List['ClinvarEvalDates']] = relationship('ClinvarEvalDates', uselist=True, back_populates='clinvar')


class DbSnp(Base):
    __tablename__ = 'db_snp'
    __table_args__ = (
        ForeignKeyConstraint(['external_db_snp_id'], ['external_references.id'], name='fk_external_references'),
        PrimaryKeyConstraint('id', name='db_snp_pkey'),
        UniqueConstraint('external_db_snp_id', name='db_snp_external_db_snp_id_key')
    )

    id = mapped_column(Integer, Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False, cache=1))
    rsid = mapped_column(String(15), nullable=False)
    external_db_snp_id = mapped_column(Integer, nullable=False)

    external_db_snp: Mapped['ExternalReferences'] = relationship('ExternalReferences', back_populates='db_snp')


class GeneAnnotations(Base):
    __tablename__ = 'gene_annotations'
    __table_args__ = (
        PrimaryKeyConstraint('id', name='gene_annotations_pkey'),
    )

    id = mapped_column(Integer,
                       Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False,
                                cache=1))
    seq_name = mapped_column(Text, nullable=False)
    source = mapped_column(Text, nullable=False)
    feature = mapped_column(Text, nullable=False)
    start_location = mapped_column(Integer, nullable=False)
    end_location = mapped_column(Integer, nullable=False)
    score = mapped_column(Double(53))
    strand = mapped_column(EnumSQL(Strand, name='strand'))
    frame = mapped_column(CHAR(1))

    gene_attributes: Mapped[List['GeneAttributes']] = relationship('GeneAttributes', uselist=True,
                                                                   back_populates='gene')
    variants: Mapped[List['Variants']] = relationship('Variants', uselist=True, back_populates='gene')


class Publications(Base):
    __tablename__ = 'publications'
    __table_args__ = (
        PrimaryKeyConstraint('id', name='publications_pkey'),
    )

    id = mapped_column(Integer,
                       Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False,
                                cache=1))
    title = mapped_column(Text)
    pmid = mapped_column(Integer)
    doi = mapped_column(Text)
    abstract = mapped_column(Text)
    match_in_sup_material = mapped_column(Boolean)
    date_published = mapped_column(Date)
    authors = mapped_column(Text)
    journal = mapped_column(Text)
    link = mapped_column(Text)

    variant: Mapped[List['Variants']] = relationship('Variants', secondary='variants_publications',
                                                     back_populates='publications')
    review: Mapped['Reviews'] = relationship('Reviews', secondary='reviews_publications', back_populates='publications')


class ScientificMembers(UserMixin, Base):
    __tablename__ = 'scientific_members'
    __table_args__ = (
        PrimaryKeyConstraint('id', name='scientific_members_pkey'),
    )

    id = mapped_column(Integer,
                       Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False,
                                cache=1))
    name = mapped_column(Text, nullable=False)
    surname = mapped_column(Text, nullable=False)
    email = mapped_column(Text, nullable=False)
    password = mapped_column(Text, nullable=False)

    variants_samples_uploads: Mapped[List['VariantsSamplesUploads']] = relationship('VariantsSamplesUploads',
                                                                                    uselist=True,
                                                                                    back_populates='scientific_member')
    reviews: Mapped[List['Reviews']] = relationship('Reviews', uselist=True, back_populates='scientific_member')


class GeneAttributes(Base):
    __tablename__ = 'gene_attributes'
    __table_args__ = (
        ForeignKeyConstraint(['gene_id'], ['gene_annotations.id'], name='fk_gene_annotations'),
        PrimaryKeyConstraint('gene_id', 'attribute_name', name='gene_attributes_pkey')
    )

    gene_id = mapped_column(Integer, nullable=False)
    attribute_name = mapped_column(Text, nullable=False)
    attribute_value = mapped_column(Text, nullable=False)

    gene: Mapped['GeneAnnotations'] = relationship('GeneAnnotations', back_populates='gene_attributes')


class Samples(Base):
    __tablename__ = 'samples'
    __table_args__ = (
        PrimaryKeyConstraint('id', name='samples_pkey'),
    )

    id = mapped_column(Text)
    genome_version = mapped_column(String(20))

    ontology_term: Mapped[List['Phenotypes']] = relationship('Phenotypes', secondary='samples_phenotypes',
                                                             back_populates='sample')
    variants_samples: Mapped[List['VariantsSamples']] = relationship('VariantsSamples', uselist=True,
                                                                     back_populates='sample')


class VariantsSamplesUploads(Base):
    __tablename__ = 'variants_samples_uploads'
    __table_args__ = (
        ForeignKeyConstraint(['variant_id', 'sample_id'], ['variants_samples.variant_id', 'variants_samples.sample_id'],
                             name='fk_variants_samples'),
        ForeignKeyConstraint(['scientific_member_id'], ['scientific_members.id'], name='fk_scientific_members'),
        PrimaryKeyConstraint('id', name='variants_samples_uploads_pkey')
    )

    id = mapped_column(Integer, Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False, cache=1))
    variant_id = mapped_column(Integer, nullable=False)
    sample_id = mapped_column(Text, nullable=False)
    upload_type = mapped_column(Text, nullable=False)
    date_uploaded = mapped_column(DateTime, nullable=False)
    scientific_member_id = mapped_column(Integer, nullable=False)

    variants_samples: Mapped['VariantsSamples'] = relationship('VariantsSamples', back_populates='variants_samples_uploads')
    scientific_member: Mapped['ScientificMembers'] = relationship('ScientificMembers', back_populates='variants_samples_uploads')
    file_upload: Mapped['FileUploads'] = relationship('FileUploads', secondary='file_uploads_variants_samples_uploads', back_populates='variants_samples_uploads')
    manual_uploads: Mapped['ManualUploads'] = relationship('ManualUploads', back_populates='variants_samples_uploads_manual')


class FileUploads(Base):
    __tablename__ = 'file_uploads'
    __table_args__ = (
        PrimaryKeyConstraint('id', name='file_uploads_pkey'),
    )

    id = mapped_column(Integer,
                       Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False,
                                cache=1))
    filename = mapped_column(Text, nullable=False)

    variants_samples_uploads: Mapped[List['VariantsSamplesUploads']] = relationship('VariantsSamplesUploads', secondary='file_uploads_variants_samples_uploads', back_populates='file_upload')


t_file_uploads_variants_samples_uploads = Table(
    'file_uploads_variants_samples_uploads', metadata,
    Column('file_upload_id', Integer, nullable=False),
    Column('variants_samples_uploads_id', Integer, nullable=False),
    ForeignKeyConstraint(['file_upload_id'], ['file_uploads.id'], name='fk_file_uploads'),
    ForeignKeyConstraint(['variants_samples_uploads_id'], ['variants_samples_uploads.id'], name='fk_variants_samples_uploads'),
    PrimaryKeyConstraint('file_upload_id', 'variants_samples_uploads_id', name='file_uploads_variants_samples_uploads_pkey')
)


class ManualUploads(Base):
    __tablename__ = 'manual_uploads'
    __table_args__ = (
        ForeignKeyConstraint(['variants_samples_uploads_manual_id'], ['variants_samples_uploads.id'], name='fk_variants_samples_uploads'),
        PrimaryKeyConstraint('id', name='manual_uploads_pkey')
    )

    id = mapped_column(Integer, Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False, cache=1))
    variants_samples_uploads_manual_id = mapped_column(Integer, nullable=False)

    variants_samples_uploads_manual: Mapped['VariantsSamplesUploads'] = relationship('VariantsSamplesUploads', back_populates='manual_uploads')


class Phenotypes(Base):
    __tablename__ = 'phenotypes'
    __table_args__ = (
        PrimaryKeyConstraint('ontology_term_id', name='phenotypes_pkey'),
    )

    ontology_term_id = mapped_column(Text)
    term_name = mapped_column(Text)

    sample: Mapped[List['Samples']] = relationship('Samples', secondary='samples_phenotypes',
                                                   back_populates='ontology_term')


t_samples_phenotypes = Table(
    'samples_phenotypes', metadata,
    Column('sample_id', Text, nullable=False),
    Column('ontology_term_id', Text, nullable=False),
    ForeignKeyConstraint(['ontology_term_id'], ['phenotypes.ontology_term_id'], name='fk_phenotypes'),
    ForeignKeyConstraint(['sample_id'], ['samples.id'], name='fk_samples'),
    PrimaryKeyConstraint('sample_id', 'ontology_term_id', name='samples_phenotypes_pkey')
)


class Variants(Base):
    __tablename__ = 'variants'
    __table_args__ = (
        ForeignKeyConstraint(['gene_id'], ['gene_annotations.id'], name='fk_gene_annotations'),
        PrimaryKeyConstraint('id', name='variants_pkey')
    )

    id = mapped_column(Integer,
                       Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False,
                                cache=1))
    chromosome = mapped_column(String(2), nullable=False)
    chromosome_position = mapped_column(Text, nullable=False)
    variant_type = mapped_column(EnumSQL(VariantType, name='variant_type'), nullable=False)
    ref = mapped_column(Text, nullable=False)
    classification = mapped_column(EnumSQL(Classification, name='classification'), nullable=False)
    gene_id = mapped_column(Integer, nullable=False)
    gene_name = mapped_column(Text, nullable=False)
    alt = mapped_column(Text)
    consequences = mapped_column(EnumSQL(Consequence, name='consequence'))

    publications: Mapped[List['Publications']] = relationship('Publications', secondary='variants_publications',
                                                              back_populates='variant')
    gene: Mapped['GeneAnnotations'] = relationship('GeneAnnotations', back_populates='variants')
    external_references: Mapped[List['ExternalReferences']] = relationship('ExternalReferences', uselist=True,
                                                                           back_populates='variant')
    reviews: Mapped[List['Reviews']] = relationship('Reviews', uselist=True, back_populates='variant')
    variants_acmg_rules: Mapped[List['VariantsAcmgRules']] = relationship('VariantsAcmgRules', uselist=True,
                                                                          back_populates='variant')
    variants_samples: Mapped[List['VariantsSamples']] = relationship('VariantsSamples', uselist=True,
                                                                     back_populates='variant')


class Reviews(Base):
    __tablename__ = 'reviews'
    __table_args__ = (
        ForeignKeyConstraint(['scientific_member_id'], ['scientific_members.id'], name='fk_scientific_members'),
        ForeignKeyConstraint(['variant_id'], ['variants.id'], name='fk_variants'),
        PrimaryKeyConstraint('id', name='reviews_pkey')
    )

    id = mapped_column(Integer,
                       Identity(always=True, start=1, increment=1, minvalue=1, maxvalue=2147483647, cycle=False,
                                cache=1))
    variant_id = mapped_column(Integer, nullable=False)
    scientific_member_id = mapped_column(Integer, nullable=False)
    date_added = mapped_column(DateTime, nullable=False)
    classification = mapped_column(EnumSQL(Classification, name='classification'), nullable=False)
    review_status = mapped_column(EnumSQL(ReviewStatus, name='review_status'))
    classification_reason = mapped_column(Text)

    publications: Mapped['Publications'] = relationship('Publications', secondary='reviews_publications',
                                                        back_populates='review')
    scientific_member: Mapped['ScientificMembers'] = relationship('ScientificMembers', back_populates='reviews')
    variant: Mapped['Variants'] = relationship('Variants', back_populates='reviews')
    acmg_rules: Mapped['AcmgRules'] = relationship('AcmgRules', secondary='reviews_acmg_rules', back_populates='review')


class VariantsAcmgRules(Base):
    __tablename__ = 'variants_acmg_rules'
    __table_args__ = (
        ForeignKeyConstraint(['acmg_rule_id'], ['acmg_rules.id'], name='fk_acmg_rules'),
        ForeignKeyConstraint(['variant_id'], ['variants.id'], name='fk_variants'),
        PrimaryKeyConstraint('variant_id', 'acmg_rule_id', name='variants_acmg_rules_pkey')
    )

    variant_id = mapped_column(Integer, nullable=False)
    acmg_rule_id = mapped_column(Integer, nullable=False)
    rule_name = mapped_column(EnumSQL(ACMGRule, name='acmg_rule'), nullable=False)

    acmg_rule: Mapped['AcmgRules'] = relationship('AcmgRules', back_populates='variants_acmg_rules')
    variant: Mapped['Variants'] = relationship('Variants', back_populates='variants_acmg_rules')


t_variants_publications = Table(
    'variants_publications', metadata,
    Column('variant_id', Integer, nullable=False),
    Column('publication_id', Integer, nullable=False),
    ForeignKeyConstraint(['publication_id'], ['publications.id'], name='fk_publications'),
    ForeignKeyConstraint(['variant_id'], ['variants.id'], name='fk_variants'),
    PrimaryKeyConstraint('variant_id', 'publication_id', name='variants_publications_pkey')
)


class VariantsSamples(Base):
    __tablename__ = 'variants_samples'
    __table_args__ = (
        ForeignKeyConstraint(['sample_id'], ['samples.id'], name='fk_samples'),
        ForeignKeyConstraint(['variant_id'], ['variants.id'], name='fk_variants'),
        PrimaryKeyConstraint('variant_id', 'sample_id', name='variants_samples_pkey')
    )

    variant_id = mapped_column(Integer, nullable=False)
    sample_id = mapped_column(Integer, nullable=False)
    genotype = mapped_column(EnumSQL(Genotype, name='genotype'), nullable=False)

    sample: Mapped['Samples'] = relationship('Samples', back_populates='variants_samples')
    variant: Mapped['Variants'] = relationship('Variants', back_populates='variants_samples')
    variants_samples_uploads: Mapped[List['VariantsSamplesUploads']] = relationship('VariantsSamplesUploads',
                                                                                    uselist=True,
                                                                                    back_populates='variants_samples')


t_reviews_acmg_rules = Table(
    'reviews_acmg_rules', metadata,
    Column('review_id', Integer, nullable=False),
    Column('acmg_rule_id', EnumSQL(ACMGRule, name='acmg_rule'), nullable=False),
    ForeignKeyConstraint(['review_id'], ['reviews.id'], name='fk_reviews'),
    ForeignKeyConstraint(['acmg_rule_id'], ['acmg_rules.id'], name='fk_acmg_rules'),
    PrimaryKeyConstraint('review_id', 'acmg_rule_id', name='reviews_acmg_rules_pkey')
)

t_reviews_publications = Table(
    'reviews_publications', metadata,
    Column('review_id', Integer, nullable=False),
    Column('publication_id', Integer, nullable=False),
    ForeignKeyConstraint(['publication_id'], ['publications.id'], name='fk_publications'),
    ForeignKeyConstraint(['review_id'], ['reviews.id'], name='fk_reviews'),
    PrimaryKeyConstraint('review_id', 'publication_id', name='reviews_publications_pkey')
)
