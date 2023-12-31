-- Create a new database
CREATE DATABASE "vus-app-db"
WITH
OWNER = postgres
ENCODING = 'UTF8'
lib = 'libc'
CONNECTION LIMIT = -1
IS_TEMPLATE = False;

-- Connect to the newly created database
\c "vus-app-db";

-- Create tables within the new database
CREATE TYPE EXTERNAL_REF_DB_TYPE AS ENUM ('DBSNP', 'CLINVAR');
CREATE TYPE VARIANT_TYPE AS ENUM ('SNV', 'MNV', 'INDEL');
CREATE TYPE CONSEQUENCE AS ENUM ('MISSENSE', 'NONSENSE', 'INSERTION', 'DELETION', 'FRAMESHIFT', 'DUPLICATION');
CREATE TYPE STRAND AS ENUM ('+', '-');
CREATE TYPE GENOTYPE AS ENUM ('HOMOZYGOUS', 'HETEROZYGOUS');
CREATE TYPE ACMG_RULE AS ENUM ('PS2', 'PM3', 'PM6', 'PP1', 'PP4', 'BS4', 'BP2', 'PS3', 'PP5', 'BP6');
CREATE TYPE ACMG_STRENGTH AS ENUM ('SUPPORTING', 'MODERATE', 'STRONG', 'VERY_STRONG');
CREATE TYPE CLASSIFICATION AS ENUM ('VUS', 'UNCERTAIN_SIGNIFICANCE', 'UNCLASSIFIED');

DROP TABLE IF EXISTS sample_files;
DROP TABLE IF EXISTS sample;
DROP TABLE IF EXISTS db_snp;
DROP TABLE IF EXISTS clinvar;
DROP TABLE IF EXISTS external_references;
DROP TABLE IF EXISTS gene_annotations;
DROP TABLE IF EXISTS variants;

-- GENE ANNOTATIONS
CREATE TABLE gene_annotations (
    gene_id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    seq_name VARCHAR(2) NOT NULL,
    source TEXT NOT NULL,
    feature TEXT NOT NULL,
    start_location INT NOT NULL,
    end_location INT NOT NULL,
    score FLOAT NOT NULL,
    strand STRAND,
    frame CHAR,
    attribute TEXT,
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
    CONSTRAINT fk_gene_annotations 
        FOREIGN KEY (gene_id) 
            REFERENCES gene_annotations(gene_id)
);

-- EXTERNAL REFERENCES
CREATE TABLE db_snp (
    db_snp_id VARCHAR(15) PRIMARY KEY
);

CREATE TABLE clinvar (
    clinvar_id TEXT PRIMARY KEY,
    canonical_spdi TEXT NOT NULL,
    classification TEXT,
    last_evaluated TIMESTAMP,
    review_status TEXT
);

CREATE TABLE external_references (
    external_references_id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    error_msg TEXT,
    db_id VARCHAR(15) NOT NULL,
    db_type EXTERNAL_REF_DB_TYPE NOT NULL,
    CONSTRAINT fk_clinvar
        FOREIGN KEY (db_id) 
            REFERENCES clinvar(clinvar_id) DEFERRABLE INITIALLY DEFERRED,
    CONSTRAINT fk_db_snp
        FOREIGN KEY (db_id) 
            REFERENCES db_snp(db_snp_id) DEFERRABLE INITIALLY DEFERRED,
    CHECK (
        (db_type = 'clinvar' AND EXISTS (SELECT 1 FROM clinvar WHERE clinvar_id = db_id)) OR
        (db_type = 'dbsnp' AND EXISTS (SELECT 1 FROM db_snp WHERE db_snp_id = db_id))
    )
);

-- VARIANTS/EXTERNAL_REFERENCES
CREATE TABLE variants_external_references ()

-- SAMPLES
CREATE TABLE sample_files (
    sample_file_id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    date_uploaded TIMESTAMP NOT NULL,
    filename TEXT NOT NULL,
    are_rsids_retrieved BOOLEAN NOT NULL,
    is_clinvar_accessed BOOLEAN NOT NULL
);

CREATE TABLE samples (
    sample_id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    phenotype TEXT,
    date_collected TIMESTAMP,
    genome_version VARCHAR(20),
    sample_file_id INT NOT NULL,
    CONSTRAINT fk_sample_files
        FOREIGN KEY (sample_file_id) 
            REFERENCES sample_files(sample_file_id)

);

-- VARIANTS/SAMPLES
CREATE TABLE variants_samples (
    variant_id INT NOT NULL,
    sample_id INT NOT NULL,
    genotype GENOTYPE NOT NULL,
    CONSTRAINT fk_variants
        FOREIGN KEY (variant_id) 
            REFERENCES variants(variant_id),
    CONSTRAINT fk_samples
        FOREIGN KEY (sample_id) 
            REFERENCES samples(sample_id)
    PRIMARY KEY (variant_id, sample_id),
) 

-- ACMG_RULES
CREATE TABLE acmg_rules(
    rule_name ACMG_RULE PRIMARY KEY,
    description TEXT NOT NULL,
    default_strength ACMG_STRENGTH NOT NULL,
    requires_lab_verification BOOLEAN NOT NULL,
)

-- SAMPLES/VARIANTS/ACMG_RULES
CREATE TABLE samples_variants_acmg_rules(
    variant_id INT NOT NULL,
    sample_id INT NOT NULL,
    rule_name ACMG_RULE NOT NULL,
    CONSTRAINT fk_variants
        FOREIGN KEY (variant_id) 
            REFERENCES variants(variant_id),
    CONSTRAINT fk_samples
        FOREIGN KEY (sample_id) 
            REFERENCES samples(sample_id)
    CONSTRAINT fk_acmg_rules
        FOREIGN KEY (rule_name) 
            REFERENCES acmg_rules(rule_name)
    PRIMARY KEY (variant_id, sample_id, rule_name),
)

-- SCIENTIFIC_MEMBERS
CREATE TABLE scientific_members(
    scientific_member_id INT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    job_title TEXT NOT NULL,
    department TEXT NOT NULL
)

-- CLASSIFICATION_OVERRIDES (VARIANTS/SCIENTIFIC_MEMBERS)
CREATE TABLE classification_overrides(
    variant_id INT NOT NULL,
    scientific_member_id INT NOT NULL,
    classification CLASSIFICATION NOT NULL,
    date_added TIMESTAMP NOT NULL,
    reason TEXT,
    CONSTRAINT fk_variants
        FOREIGN KEY (variant_id) 
            REFERENCES variants(variant_id),
    CONSTRAINT fk_scientific_members
        FOREIGN KEY (scientific_member_id) 
            REFERENCES scientific_members(scientific_member_id)
    PRIMARY KEY (variant_id, scientific_member_id),
)

-- CLASSIFICATION_OVERRIDES/ACMG_RULES
CREATE TABLE classification_overrides_acmg_rules(
    variant_id INT NOT NULL,
    scientific_member_id INT NOT NULL,
    rule_name ACMG_RULE NOT NULL,
    CONSTRAINT fk_variants
        FOREIGN KEY (variant_id) 
            REFERENCES variants(variant_id),
    CONSTRAINT fk_scientific_members
        FOREIGN KEY (scientific_member_id) 
            REFERENCES scientific_members(scientific_member_id)
    CONSTRAINT fk_acmg_rules
        FOREIGN KEY (rule_name) 
            REFERENCES acmg_rules(rule_name)
    PRIMARY KEY (variant_id, scientific_member_id, rule_name),
)

-- PUBLICATIONS
CREATE TABLE publications(
    pmid INT PRIMARY KEY,
    title TEXT NOT NULL,
    abstract TEXT,
    match_in_sup_material BOOLEAN NOT NULL,
    date_published DATE NOT NULL
)

-- VARIANTS/PUBLICATIONS
CREATE TABLE variants_publications(
    variant_id INT NOT NULL,
    pmid INT NOT NULL,
    CONSTRAINT fk_variants
        FOREIGN KEY (variant_id) 
            REFERENCES variants(variant_id),
    CONSTRAINT fk_publications
        FOREIGN KEY (pmid) 
            REFERENCES publications(pmid)
    PRIMARY KEY (variant_id, pmid),
)

-- CLASSIFICATION_OVERRIDES/PUBLICATIONS
CREATE TABLE classification_overrides_publications(
    variant_id INT NOT NULL,
    scientific_member_id INT NOT NULL,
    pmid INT NOT NULL,
    CONSTRAINT fk_variants
        FOREIGN KEY (variant_id) 
            REFERENCES variants(variant_id),
    CONSTRAINT fk_scientific_members
        FOREIGN KEY (scientific_member_id) 
            REFERENCES scientific_members(scientific_member_id)
    CONSTRAINT fk_publications
        FOREIGN KEY (pmid) 
            REFERENCES publications(pmid)
    PRIMARY KEY (variant_id, scientific_member_id, pmid),
)




