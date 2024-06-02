-- First run this:
-- Create a new database
CREATE DATABASE "vus-app-db"
WITH
OWNER = postgres
ENCODING = 'UTF8'
CONNECTION LIMIT = -1
IS_TEMPLATE = False;


-- Once db is created, change connection to new db and run this:
-- Create tables within the new database
-- CREATE TYPE EXTERNAL_REF_DB_TYPE AS ENUM ('DBSNP', 'CLINVAR');
CREATE TYPE VARIANT_TYPE AS ENUM ('SNV', 'MNV', 'INDEL');
CREATE TYPE CONSEQUENCE AS ENUM ('MISSENSE', 'NONSENSE', 'INSERTION', 'DELETION', 'FRAMESHIFT', 'DUPLICATION');
CREATE TYPE STRAND AS ENUM ('POSITIVE', 'NEGATIVE');
CREATE TYPE GENOTYPE AS ENUM ('HOMOZYGOUS', 'HETEROZYGOUS');
CREATE TYPE ACMG_RULE AS ENUM ('PS2', 'PM3', 'PM6', 'PP1', 'PP4', 'BS4', 'BP2', 'PS3', 'PP5', 'BP6', 'PVS1', 'PS1', 'BS3', 'PM1', 'BP3', 'PM2', 'PM4', 'PM5', 'PP2', 'BP1', 'PP3', 'BP4', 'BA1', 'BS1', 'BS2', 'BP7', 'PS4', 'BP5');
CREATE TYPE ACMG_STRENGTH AS ENUM ('SUPPORTING', 'MODERATE', 'STRONG', 'VERY_STRONG', 'STAND_ALONE');
CREATE TYPE CLASSIFICATION AS ENUM ('PATHOGENIC', 'LIKELY_PATHOGENIC', 'VUS', 'LIKELY_BENIGN', 'BENIGN');
CREATE TYPE REVIEW_STATUS AS ENUM ('IN_PROGRESS', 'COMPLETE');


-- DROP TABLE IF EXISTS reviews_publications;
-- DROP TABLE IF EXISTS variants_publications;
-- DROP TABLE IF EXISTS publications;
-- DROP TABLE IF EXISTS reviews_acmg_rules;
-- DROP TABLE IF EXISTS classification_overrides;
-- DROP TABLE IF EXISTS reviews;
-- DROP TABLE IF EXISTS scientific_members;
-- DROP TABLE IF EXISTS samples_variants_acmg_rules;
-- DROP TABLE IF EXISTS acmg_rules;
-- DROP TABLE IF EXISTS variants_samples;
-- DROP TABLE IF EXISTS samples;
-- DROP TABLE IF EXISTS sample_files;
-- DROP TABLE IF EXISTS variants_external_references;
-- DROP TABLE IF EXISTS external_references;
-- DROP TABLE IF EXISTS clinvar;
-- DROP TABLE IF EXISTS db_snp;
-- DROP TABLE IF EXISTS variants;
-- DROP TABLE IF EXISTS gene_annotations;


-- GENE ANNOTATIONS
CREATE TABLE gene_annotations (
    id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    seq_name TEXT NOT NULL,
    source TEXT NOT NULL,
    feature TEXT NOT NULL,
    start_location INT NOT NULL,
    end_location INT NOT NULL,
    score FLOAT,
    strand STRAND,
    frame CHAR
);

CREATE TABLE gene_attributes (
    gene_id INT NOT NULL,
    attribute_name TEXT NOT NULL,
    attribute_value TEXT NOT NULL,
    CONSTRAINT fk_gene_annotations
        FOREIGN KEY (gene_id) 
            REFERENCES gene_annotations(id),
    PRIMARY KEY (gene_id, attribute_name)
);

-- VARIANTS 
CREATE TABLE variants (
    id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    chromosome VARCHAR(2) NOT NULL,
    chromosome_position TEXT NOT NULL,
    variant_type VARIANT_TYPE NOT NULL,
    ref TEXT NOT NULL,
    alt TEXT,
    consequences CONSEQUENCE,
    classification CLASSIFICATION NOT NULL,
    gene_id INT NOT NULL,
    gene_name TEXT NOT NULL,
    CONSTRAINT fk_gene_annotations 
        FOREIGN KEY (gene_id) 
            REFERENCES gene_annotations(id)
);

-- EXTERNAL REFERENCES
-- Table that references either clinvar or db_snp based on db_type
CREATE TABLE external_references (
    id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    error_msg TEXT,
    db_type TEXT NOT NULL,
    variant_id INT NOT NULL,
    CONSTRAINT fk_variants
            FOREIGN KEY (variant_id) 
                REFERENCES variants(id)
);

CREATE TABLE db_snp (
    id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    rsid VARCHAR(15) NOT NULL,
    external_db_snp_id INT NOT NULL UNIQUE,
    CONSTRAINT fk_external_references
            FOREIGN KEY (external_db_snp_id) 
                REFERENCES external_references(id)
);

CREATE TABLE clinvar (
    id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    variation_id TEXT NOT NULL,
    external_clinvar_id INT NOT NULL UNIQUE,
    canonical_spdi TEXT NOT NULL,
    CONSTRAINT fk_external_references
        FOREIGN KEY (external_clinvar_id) 
            REFERENCES external_references(id)
);

CREATE TABLE auto_clinvar_updates (
    id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    classification TEXT,
    review_status TEXT,
	last_evaluated TIMESTAMP
);

CREATE TABLE auto_clinvar_eval_dates (
    id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    clinvar_id INT NOT NULL,
	auto_clinvar_update_id INT,
    eval_date TIMESTAMP,
    CONSTRAINT fk_clinvar
        FOREIGN KEY (clinvar_id) 
            REFERENCES clinvar(id),
	CONSTRAINT fk_auto_clinvar_updates
        FOREIGN KEY (auto_clinvar_update_id) 
            REFERENCES auto_clinvar_updates(id)
);

-- SCIENTIFIC_MEMBERS
CREATE TABLE scientific_members(
    id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    name TEXT NOT NULL,
    surname TEXT NOT NULL,
	email TEXT NOT NULL,
	password TEXT NOT NULL
);

-- VARIANT HGVS
CREATE TABLE variant_hgvs(
    id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
	variant_id INT NOT NULL,
    hgvs TEXT NOT NULL,
    CONSTRAINT fk_variants
        FOREIGN KEY (variant_id) 
            REFERENCES variants(id)
);

-- SAMPLES
CREATE TABLE samples (
    id TEXT PRIMARY KEY,
    -- date_collected TIMESTAMP,
    genome_version VARCHAR(20)
);

CREATE TABLE phenotypes (
  ontology_term_id TEXT PRIMARY KEY,
  term_name TEXT  
);

-- SAMPLES/PHENOTYPES
CREATE TABLE samples_phenotypes(
    sample_id TEXT NOT NULL,
    ontology_term_id TEXT NOT NULL,
    CONSTRAINT fk_samples
        FOREIGN KEY (sample_id) 
            REFERENCES samples(id),
    CONSTRAINT fk_phenotypes
        FOREIGN KEY (ontology_term_id) 
            REFERENCES phenotypes(ontology_term_id),
    PRIMARY KEY (sample_id, ontology_term_id)
);

-- VARIANTS/SAMPLES
CREATE TABLE variants_samples (
    variant_id INT NOT NULL,
    sample_id TEXT NOT NULL,
	variant_hgvs_id INT NOT NULL,
    genotype GENOTYPE NOT NULL,
	CONSTRAINT fk_variant_hgvs
        FOREIGN KEY (variant_hgvs_id) 
            REFERENCES variant_hgvs(id),
    CONSTRAINT fk_variants
        FOREIGN KEY (variant_id) 
            REFERENCES variants(id),
    CONSTRAINT fk_samples
        FOREIGN KEY (sample_id) 
            REFERENCES samples(id),
    PRIMARY KEY (variant_id, sample_id)
); 


-- Table that references either file or manual upload based on upload_type
-- VARIANTS SAMPLES UPLOADS
CREATE TABLE variants_samples_uploads (
    id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
	variant_id INT NOT NULL,
    sample_id TEXT NOT NULL,
    upload_type TEXT NOT NULL,
	date_uploaded TIMESTAMP NOT NULL,
    scientific_member_id INT NOT NULL,
     CONSTRAINT fk_variants_samples
        FOREIGN KEY (variant_id, sample_id) 
            REFERENCES variants_samples(variant_id, sample_id),
    CONSTRAINT fk_scientific_members
        FOREIGN KEY (scientific_member_id) 
            REFERENCES scientific_members(id)
);

CREATE TABLE file_uploads (
    id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    filename TEXT NOT NULL
    -- are_rsids_retrieved BOOLEAN NOT NULL,
    -- is_clinvar_accessed BOOLEAN NOT NULL
);


-- FILE UPLOADS/VARIANTS SAMPLES UPLOADS
CREATE TABLE file_uploads_variants_samples_uploads(
    file_upload_id INT NOT NULL,
    variants_samples_uploads_id INT NOT NULL,
    CONSTRAINT fk_file_uploads
        FOREIGN KEY (file_upload_id) 
            REFERENCES file_uploads(id),
    CONSTRAINT fk_variants_samples_uploads
        FOREIGN KEY (variants_samples_uploads_id) 
            REFERENCES variants_samples_uploads(id),
    PRIMARY KEY (file_upload_id, variants_samples_uploads_id)
);


CREATE TABLE manual_uploads (
    id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    variants_samples_uploads_manual_id INT NOT NULL,
    CONSTRAINT fk_variants_samples_uploads
            FOREIGN KEY (variants_samples_uploads_manual_id) 
                REFERENCES variants_samples_uploads(id)
);

-- ACMG_RULES
CREATE TABLE acmg_rules(
    id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    rule_name ACMG_RULE,
    description TEXT NOT NULL,
    default_strength ACMG_STRENGTH NOT NULL,
    requires_lab_verification BOOLEAN NOT NULL
);

-- VARIANTS/ACMG_RULES
CREATE TABLE variants_acmg_rules(
    variant_id INT NOT NULL,
    acmg_rule_id INT NOT NULL,
    rule_name ACMG_RULE NOT NULL,
    CONSTRAINT fk_variants
        FOREIGN KEY (variant_id) 
            REFERENCES variants(id),
    CONSTRAINT fk_acmg_rules
        FOREIGN KEY (acmg_rule_id) 
            REFERENCES acmg_rules(id),
    PRIMARY KEY (variant_id, acmg_rule_id)
);

-- REVIEWS (VARIANTS/SCIENTIFIC_MEMBERS)
CREATE TABLE reviews(
    id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    variant_id INT NOT NULL,
    scientific_member_id INT NOT NULL,
    date_added TIMESTAMP NOT NULL,
    review_status REVIEW_STATUS,
    classification CLASSIFICATION NOT NULL,
    classification_reason TEXT,
	is_acmg_rule_added BOOL,
	is_acmg_rule_deleted BOOL,
    CONSTRAINT fk_variants
        FOREIGN KEY (variant_id) 
            REFERENCES variants(id),
    CONSTRAINT fk_scientific_members
        FOREIGN KEY (scientific_member_id) 
            REFERENCES scientific_members(id)
);

-- REVIEWS/ACMG_RULES
CREATE TABLE reviews_acmg_rules(
    review_id INT NOT NULL,
    acmg_rule_id INT NOT NULL,
    CONSTRAINT fk_reviews
        FOREIGN KEY (review_id) 
            REFERENCES reviews(id),
    CONSTRAINT fk_acmg_rules
        FOREIGN KEY (acmg_rule_id) 
            REFERENCES acmg_rules(id),
    PRIMARY KEY (review_id, acmg_rule_id)
);

-- PUBLICATIONS
CREATE TABLE publications(
    id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY, 
    pmid INT,
    doi TEXT,
    title TEXT,
    abstract TEXT,
    match_in_sup_material BOOLEAN,
    date_published DATE,
	authors TEXT,
	journal TEXT,
    link TEXT
);

-- VARIANTS/PUBLICATIONS
CREATE TABLE variants_publications(
    variant_id INT NOT NULL,
    publication_id INT NOT NULL,
    CONSTRAINT fk_variants
        FOREIGN KEY (variant_id) 
            REFERENCES variants(id),
    CONSTRAINT fk_publications
        FOREIGN KEY (publication_id) 
            REFERENCES publications(id),
    PRIMARY KEY (variant_id, publication_id)
);

-- REVIEWS/PUBLICATIONS
CREATE TABLE reviews_publications(
    review_id INT NOT NULL,
    publication_id INT NOT NULL,
    CONSTRAINT fk_reviews
        FOREIGN KEY (review_id) 
            REFERENCES reviews(id),
    CONSTRAINT fk_publications
        FOREIGN KEY (publication_id) 
            REFERENCES publications(id),
    PRIMARY KEY (review_id, publication_id)
);




