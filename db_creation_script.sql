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
CREATE TYPE ACMG_RULE AS ENUM ('PS2', 'PM3', 'PM6', 'PP1', 'PP4', 'BS4', 'BP2', 'PS3', 'PP5', 'BP6');
CREATE TYPE ACMG_STRENGTH AS ENUM ('SUPPORTING', 'MODERATE', 'STRONG', 'VERY_STRONG');
CREATE TYPE CLASSIFICATION AS ENUM ('PATHOGENIC', 'LIKELY_PATHOGENIC', 'VUS', 'UNCERTAIN_SIGNIFICANCE', 'UNCLASSIFIED'. 'LIKELY_VUS', 'LIKELY_BENIGN', 'BENIGN');
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
    gene_id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
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
            REFERENCES gene_annotations(gene_id),
    PRIMARY KEY (gene_id, attribute_name)
);

-- VARIANTS 
CREATE TABLE variants (
    variant_id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
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
            REFERENCES gene_annotations(gene_id)
);

-- EXTERNAL REFERENCES
-- Table that references either clinvar or db_snp based on db_type
CREATE TABLE external_references (
    external_references_id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    error_msg TEXT,
    db_type TEXT NOT NULL,
    variant_id INT NOT NULL,
    CONSTRAINT fk_variants
            FOREIGN KEY (variant_id) 
                REFERENCES variants(variant_id)
);

CREATE TABLE db_snp (
    db_snp_id VARCHAR(15) PRIMARY KEY,
    external_db_snp_id INT NOT NULL UNIQUE,
    CONSTRAINT fk_external_references
            FOREIGN KEY (external_db_snp_id) 
                REFERENCES external_references(external_references_id)
);

CREATE TABLE clinvar (
    clinvar_id TEXT PRIMARY KEY,
    external_clinvar_id INT NOT NULL UNIQUE,
    canonical_spdi TEXT NOT NULL,
    classification TEXT,
    last_evaluated TIMESTAMP,
    review_status TEXT,
    CONSTRAINT fk_external_references
        FOREIGN KEY (external_clinvar_id) 
            REFERENCES external_references(external_references_id)
);

-- SAMPLES
CREATE TABLE sample_files (
    sample_file_id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    date_uploaded TIMESTAMP NOT NULL,
    filename TEXT NOT NULL
    -- are_rsids_retrieved BOOLEAN NOT NULL,
    -- is_clinvar_accessed BOOLEAN NOT NULL
);

CREATE TABLE samples (
    sample_id TEXT PRIMARY KEY,
    phenotype TEXT,
    -- date_collected TIMESTAMP,
    genome_version VARCHAR(20),
    sample_file_id INT NOT NULL,
    CONSTRAINT fk_sample_files
        FOREIGN KEY (sample_file_id) 
            REFERENCES sample_files(sample_file_id)

);

-- VARIANTS/SAMPLES
CREATE TABLE variants_samples (
    variant_id INT NOT NULL,
    sample_id TEXT NOT NULL,
    genotype GENOTYPE NOT NULL,
    CONSTRAINT fk_variants
        FOREIGN KEY (variant_id) 
            REFERENCES variants(variant_id),
    CONSTRAINT fk_samples
        FOREIGN KEY (sample_id) 
            REFERENCES samples(sample_id),
    PRIMARY KEY (variant_id, sample_id)
); 

-- ACMG_RULES
CREATE TABLE acmg_rules(
    rule_name ACMG_RULE PRIMARY KEY,
    description TEXT NOT NULL,
    default_strength ACMG_STRENGTH NOT NULL,
    requires_lab_verification BOOLEAN NOT NULL
);

-- SAMPLES/VARIANTS/ACMG_RULES
CREATE TABLE samples_variants_acmg_rules(
    variant_id INT NOT NULL,
    sample_id TEXT NOT NULL,
    rule_name ACMG_RULE NOT NULL,
    CONSTRAINT fk_variants
        FOREIGN KEY (variant_id) 
            REFERENCES variants(variant_id),
    CONSTRAINT fk_samples
        FOREIGN KEY (sample_id) 
            REFERENCES samples(sample_id),
    CONSTRAINT fk_acmg_rules
        FOREIGN KEY (rule_name) 
            REFERENCES acmg_rules(rule_name),
    PRIMARY KEY (variant_id, sample_id, rule_name)
);

-- SCIENTIFIC_MEMBERS
CREATE TABLE scientific_members(
    scientific_member_id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    name TEXT NOT NULL,
    surname TEXT NOT NULL
);

-- REVIEWS (VARIANTS/SCIENTIFIC_MEMBERS)
CREATE TABLE reviews(
    review_id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    variant_id INT NOT NULL,
    scientific_member_id INT NOT NULL,
    date_added TIMESTAMP NOT NULL,
    review_status REVIEW_STATUS,
    classification CLASSIFICATION NOT NULL,
    classification_reason TEXT,
    CONSTRAINT fk_variants
        FOREIGN KEY (variant_id) 
            REFERENCES variants(variant_id),
    CONSTRAINT fk_scientific_members
        FOREIGN KEY (scientific_member_id) 
            REFERENCES scientific_members(scientific_member_id)
);

-- REVIEWS/ACMG_RULES
CREATE TABLE reviews_acmg_rules(
    review_id INT NOT NULL,
    rule_name ACMG_RULE NOT NULL,
    CONSTRAINT fk_reviews
        FOREIGN KEY (review_id) 
            REFERENCES reviews(review_id),
    CONSTRAINT fk_acmg_rules
        FOREIGN KEY (rule_name) 
            REFERENCES acmg_rules(rule_name),
    PRIMARY KEY (review_id, rule_name)
);

-- PUBLICATIONS
CREATE TABLE publications(
    pmid INT PRIMARY KEY,
    title TEXT NOT NULL,
    abstract TEXT,
    match_in_sup_material BOOLEAN NOT NULL,
    date_published DATE NOT NULL
);

-- VARIANTS/PUBLICATIONS
CREATE TABLE variants_publications(
    variant_id INT NOT NULL,
    pmid INT NOT NULL,
    CONSTRAINT fk_variants
        FOREIGN KEY (variant_id) 
            REFERENCES variants(variant_id),
    CONSTRAINT fk_publications
        FOREIGN KEY (pmid) 
            REFERENCES publications(pmid),
    PRIMARY KEY (variant_id, pmid)
);

-- REVIEWS/PUBLICATIONS
CREATE TABLE reviews_publications(
    review_id INT NOT NULL,
    pmid INT NOT NULL,
    CONSTRAINT fk_reviews
        FOREIGN KEY (review_id) 
            REFERENCES reviews(review_id),
    CONSTRAINT fk_publications
        FOREIGN KEY (pmid) 
            REFERENCES publications(pmid),
    PRIMARY KEY (review_id, pmid)
);




