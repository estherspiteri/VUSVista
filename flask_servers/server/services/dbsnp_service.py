import re

from flask import current_app, Response
import pandas as pd
from Bio import Entrez
import requests
from requests import RequestException

from server.responses.internal_response import InternalResponse


def get_nucleotide_seq(chr: int, start_position: int, end_position: int):
    current_app.logger.error(
        f'Retrieving nucleotide sequence at Chr {chr} from position {start_position} to {end_position}')
    url = f"https://rest.ensembl.org/sequence/region/human/{chr}:{start_position}..{end_position}:1?content-type=text/plain&coord_system_version=GRCh37"

    try:
        get_nucleotide_sequence_res = requests.get(url)
    except RequestException as e:
        current_app.logger.error(f'Failed to connect to Ensembl API: {e}')
        # send service unavailable status code
        return InternalResponse(None, 503, e)

    if get_nucleotide_sequence_res.status_code != 200:
        current_app.logger.error(f'Ensembl API failed: {get_nucleotide_sequence_res.reason}')
        return InternalResponse(None, get_nucleotide_sequence_res.status_code, get_nucleotide_sequence_res.reason)
    else:
        seq = get_nucleotide_sequence_res.text

        return InternalResponse(seq, 200)


# convert variant list to VCF format
def convert_variants_to_vcf(variant_df: pd.DataFrame, variants_vcf_filename: str):
    with open(variants_vcf_filename, 'w') as vcf_f:
        for var in variant_df.iterrows():
            # for insertions, if no ref allele is provided, use nucleotide at given position
            if 'insertion' in var[1]['Type'].lower() and var[1]['Reference'] is None:
                get_nucleotide_seq_res = get_nucleotide_seq(var[1]['Chr'], var[1]['Position'], var[1]['Position'])

                ref = '.'
                alt = var[1]['Alt']
                if get_nucleotide_seq_res.status == 200:
                    ref = get_nucleotide_seq_res.data
                    alt = get_nucleotide_seq_res.data + alt

                vcf_string = f"{var[1]['Chr']} {var[1]['Position']} . {ref} {alt} . PASS\n"
            # for deletions, if no alt allele is provided, use nucleotide just before given position
            # TODO: cater for position 0
            elif 'deletion' in var[1]['Type'].lower() and var[1]['Alt'] is None and int(var[1]['Position']) != 0:
                get_nucleotide_seq_res = get_nucleotide_seq(var[1]['Chr'], int(var[1]['Position']) - 1,
                                                            int(var[1]['Position']) - 1)

                ref = var[1]['Reference']
                alt = '.'
                position = var[1]['Position']
                if get_nucleotide_seq_res.status == 200:
                    ref = get_nucleotide_seq_res.data + ref
                    alt = get_nucleotide_seq_res.data
                    position = int(position) - 1

                vcf_string = f"{var[1]['Chr']} {position} . {ref} {alt} . PASS\n"
            else:
                vcf_string = f"{var[1]['Chr']} {var[1]['Position']} . {var[1]['Reference']} {var[1]['Alt']} . PASS\n"

            vcf_f.write(vcf_string)


def get_rsids(genome_version: str, variants_vcf_filename: str) -> InternalResponse:
    # retrieve RSIDs if they exist for a given variant
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/vcf/file/set_rsids?assembly={genome_version}"
    myfiles = {'file': open(variants_vcf_filename, 'rb')}

    try:
        rsid_vcf_res = requests.post(url, files=myfiles)
    except RequestException as e:
        current_app.logger.error(f'Failed to connect to NCBI Variation Service: {e}')
        # send service unavailable status code
        return InternalResponse(None, 503, e)

    if rsid_vcf_res.status_code != 200:
        current_app.logger.error(f'NCBI Variation Service failed: {rsid_vcf_res.reason}')
        return InternalResponse(None, rsid_vcf_res.status_code, rsid_vcf_res.reason)
    else:
        # TODO: save updated file with rsids (overwrite initially created file)?
        rsids = []

        for variant_rsid_res in rsid_vcf_res.text.split('PASS\n'):
            if re.match("^# Error in the next line: .*\n\d+\t\d+\t\.\t\w\t\w\t\.", variant_rsid_res,
                        re.IGNORECASE | re.DOTALL):
                # TODO: try and get RSID with refAllele set to 'N'
                rsids.append('NORSID')
            # skip empty lines (x.strip()) - trying to split an empty line can lead to an "index out of range" error
            elif variant_rsid_res.strip():
                # extract RSIDs
                rsids.append(
                    variant_rsid_res.split()[2])  # columns: 'chromosome', 'position', 'rsid', 'reference', 'alt'

        return InternalResponse(rsids, 200)


# Function to get variant info from dbSNP
def get_dbsnp_variant_info(rsid: str) -> InternalResponse:
    variant_id = None
    variant_record = None

    try:
        # look up rsid in dbSNP
        search_results_handle = Entrez.esearch(db="snp", term=rsid)
        search_results = Entrez.read(search_results_handle)
        search_results_handle.close()

        # assuming top most result is most relevant
        if len(search_results['IdList']) > 0:
            variant_id = search_results['IdList'][0]
    except IOError as e:
        current_app.logger.error(f'Network error when calling Entrez.esearch(): {str(e)}')
        return InternalResponse(None, e.errno, str(e))

    if variant_id is not None:
        try:
            # retrieving most relevant variant's info
            variant_handle = Entrez.esummary(db="snp", id=variant_id)
            variant_record = Entrez.read(variant_handle)
            variant_handle.close()
        except IOError as e:
            current_app.logger.error(f'Network error when calling Entrez.efetch(): {str(e)}')
            return InternalResponse(None, e.errno, str(e))

    return InternalResponse(variant_record, 200)


# Function to get genes from dbSNP variant info
def get_genes_from_dbsnp_info(dbsnp_info):
    gene_names = []

    if dbsnp_info:
        # assuming first document summary is the most relevant
        genes = dbsnp_info['DocumentSummarySet']['DocumentSummary'][0]['GENES']

        gene_names = [gene['NAME'] for gene in genes]

    return gene_names


# Function to get chromosome and chromosome position from dbSNP variant info
def get_chr_pos_from_dbsnp_info(dbsnp_info):
    chr = ''
    pos = ''

    if dbsnp_info:
        # assuming first document summary is the most relevant
        chr_pos = dbsnp_info['DocumentSummarySet']['DocumentSummary'][0]['CHRPOS_PREV_ASSM']
        chr_pos_split = chr_pos.split(':')

        chr = chr_pos_split[0]
        pos = chr_pos_split[1]

    return {'CHR': chr, 'POS': pos}


# Function to get Reference and Alt from dbSNP - ONLY WORKS FOR SNV aka REF>ALT
def get_alleles_from_dbsnp_info(dbsnp_info, pos, variant_type: str):
    ref_alt = []
    doc_sum_array_cleaned = []
    matching_pos_doc_sum = []

    if dbsnp_info:
        # assuming first document summary is the most relevant
        doc_sum_array = dbsnp_info['DocumentSummarySet']['DocumentSummary'][0]['DOCSUM'].split(',')

        for doc_sum in doc_sum_array:
            hgvs = doc_sum

            # TODO: use regex
            if '=' in hgvs:
                hgvs = hgvs.split('=')[1]

            if pos in hgvs:
                matching_pos_doc_sum.append(hgvs)

            doc_sum_array_cleaned.append(hgvs)

        for doc_sum in matching_pos_doc_sum:
            # parsing hgvs notation to access reference and alt
            ref_alt_allele = ''.join([x for x in doc_sum.split(':')[1].split('.')[1] if not x.isdigit()])

            if variant_type == 'SNV':
                if '>' in ref_alt_allele:
                    ref_alt_allele_split = ref_alt_allele.split('>')

                    ref_alt.append({'REF': ref_alt_allele_split[0], 'ALT': ref_alt_allele_split[1]})
            # elif 'insert' in variant_type.lower():
            # elif 'deletion' in variant_type.lower():

    return ref_alt, doc_sum_array_cleaned


# Function to verify that the rsid matches with the requested variant using dbSNP
def verify_rsid(rsid: str, genes: str, chr: str, pos: str, ref: str, alt: str, type: str, user_rsid: str) -> InternalResponse:
    is_valid = True
    error_msgs = []

    db_snp_variant_info_res = get_dbsnp_variant_info(rsid)

    if db_snp_variant_info_res.status != 200:
        error_msg = f'Retrieval of variant info from dbSNP failed!'
        current_app.logger.error(error_msg)
        error_msgs.append(error_msg)
        return InternalResponse(None, 500)
    else:
        db_snp_variant_info = db_snp_variant_info_res.data

        if db_snp_variant_info is not None:
            # compare with any provided RSID
            if str(user_rsid) != "nan" and len(user_rsid) > 0 and user_rsid != rsid:
                error_msg = f"Variant's provided RSID {user_rsid} does not match the suggested RSID {rsid}!"
                current_app.logger.warn(error_msg)
                error_msgs.append(error_msg)
                is_valid = False

            # compare genes
            genes_list = get_genes_from_dbsnp_info(db_snp_variant_info)

            # flag that indicates whether at least one of the variant's genes matched the RSID's genes
            is_gene_found = False
            for gene in genes.split(','):
                if gene in genes_list:
                    is_gene_found = True
                    break

            if not is_gene_found:
                error_msg = f"Variant's genes {genes} do not match the RSID's genes {genes_list}!"
                current_app.logger.warn(error_msg)
                error_msgs.append(error_msg)
                is_valid = False

            # compare chromosome and chromosome position
            chr_pos = get_chr_pos_from_dbsnp_info(db_snp_variant_info)

            if chr_pos['CHR'] != chr:
                error_msg = f"Variant's chromosome {chr} does not match RSID's chromosome {chr_pos['CHR']}!"
                current_app.logger.warn(error_msg)
                error_msgs.append(error_msg)
                is_valid = False

            if chr_pos['POS'] != pos:
                error_msg = f"Variant's chromosome position {pos} does not match RSID's chromosome position {chr_pos['POS']}!"
                current_app.logger.warn(error_msg)
                error_msgs.append(error_msg)
                is_valid = False

            # compare reference and Alt
            ref_alt, doc_sum_array_cleaned = get_alleles_from_dbsnp_info(db_snp_variant_info, chr_pos['POS'], type)

            ref_alt_allele_match = False

            for ref_alt_allele in ref_alt:
                if ref_alt_allele['REF'] == ref and ref_alt_allele['ALT'] == alt:
                    ref_alt_allele_match = True
                    break

            if not ref_alt_allele_match:
                error_msg = (f"Variant's reference {ref} and/or alt allele {alt} do not match any of the reference and "
                             f"alt alleles of the RSID's associated HGVS")  #: {doc_sum_array_cleaned}!"
                current_app.logger.warn(error_msg)
                error_msgs.append(error_msg)
                is_valid = False
        else:
            error_msg = f"No variant info found in dbSNP!"
            current_app.logger.warn(error_msg)
            error_msgs.append(error_msg)
            is_valid = False

        return InternalResponse({'isValid': is_valid, 'errorMsgs': error_msgs}, 200)


def get_rsids_from_dbsnp(vus_df: pd.DataFrame) -> InternalResponse:
    # generate VCF string for VUS
    convert_variants_to_vcf(vus_df, 'variants.vcf')

    get_rsids_res: InternalResponse = get_rsids('GRCh37.p13', 'variants.vcf')

    if get_rsids_res.status != 200:
        current_app.logger.error(
            f'Get RSIDs query failed {get_rsids_res.status}')
        return InternalResponse(None, 500)
    else:
        # get variant RSIDs
        vus_df['RSID'] = get_rsids_res.data

        # check validity of RSIDs
        rsid_verification = []

        for index, row in vus_df.iterrows():
            if row['RSID'] == 'NORSID':
                rsid_verification.append({'isValid': False, 'errorMsgs': []})
            else:
                try:
                    verify_rsid_res = verify_rsid(row['RSID'], row['Gene'], row['Chr'], row['Position'], row['Reference'],
                                                  row['Alt'], row['Type'], row['RSID_'])
                except Exception as e:
                    print('here', e)

                if verify_rsid_res.status != 200:
                    current_app.logger.error(f"RSID verification for {row['RSID']}  failed 500")
                    return InternalResponse(None, 500)
                else:
                    rsid_verification.append(verify_rsid_res.data)

        # create a column which shows if an rsid is verified successfully or not (rsid variant matches details of inputted variant)
        vus_df['RSID dbSNP verified'] = [x['isValid'] for x in rsid_verification]
        # create a column which shows the error messages that where returned during RSID verification§
        vus_df['RSID dbSNP errorMsgs'] = ['|| '.join(x['errorMsgs']) for x in rsid_verification]

        return InternalResponse(vus_df, 200)
