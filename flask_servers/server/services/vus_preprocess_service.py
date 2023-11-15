from flask import current_app, Response
import pandas as pd
from Bio import Entrez
from werkzeug.datastructures import FileStorage
import json

from server.responses.internal_response import InternalResponse
from server.services.dbsnp_service import get_rsids_from_dbsnp

from server.services.clinvar_service import retrieve_clinvar_variant_classifications

Entrez.email = "esther.spiteri.18@um.edu.mt"


# generate separate records for each gene when a VUS has multiple genes
def handle_vus_with_multiple_genes(vus_df: pd.DataFrame):
    multiple_genes = []
    rows_to_delete = []
    new_rows_dictionaries = []
    var_id = 0

    # insert id column so that variants with multiple genes have multiple entries in the dataframe with the same id
    vus_df.insert(0, 'VUS Id', 0)

    # make a copy of the dataframe to be able to iterate through it whilst modifying the original dataframe
    new_vus_df = vus_df.copy()

    # note VUS with multiple genes and generate records for each gene (same VUS ID)
    for index, row in new_vus_df.iterrows():
        # set variant id
        vus_df.at[index, 'VUS Id'] = var_id

        # get the gene name
        gene_string = row['Gene']
        genes = gene_string.split(',')

        # note if variant has multiple genes
        filtered_genes = [g for g in genes if 'AS1' not in g and 'LOC' not in g]
        if len(filtered_genes) > 1:
            multiple_genes.append(var_id)

        if len(genes) > 1:
            rows_to_delete.append(index)

            # create a variant record for every gene
            for i, gene in enumerate(genes):
                # exclude anti-sense and locus
                if 'AS1' not in gene and 'LOC' not in gene:
                    new_row = row.copy().to_dict()
                    new_row['Gene'] = gene
                    new_row['VUS Id'] = var_id
                    new_rows_dictionaries.append(new_row)

        # increment variant id
        var_id += 1

        # TODO: handle cases where not SNVs

    current_app.logger.info(f"Percentage of VUS with multiple genes: {round(len(multiple_genes) / var_id * 100, 2)}%\n"
                            f"VUS IDs of VUS with multiple genes: {[g for g in multiple_genes]}")

    # clean up - remove multiple gene records (since they will be replaced by individual records)
    # drop rows containing variants with multiple genes (including those with AS1 and LOC)
    vus_df = vus_df.drop([g for g in rows_to_delete])
    vus_df = vus_df.reset_index(drop=True)

    # insert row for every one of the multiple genes
    for r in new_rows_dictionaries:
        df_dictionary = pd.DataFrame([r])
        vus_df = pd.concat([vus_df, df_dictionary], ignore_index=True)

    # sort the dataframe by the 'VUS Id' column in ascending order
    vus_df = vus_df.sort_values(by='VUS Id')

    # reset the indices to maintain a continuous index
    vus_df = vus_df.reset_index(drop=True)

    return vus_df


# create columns for chromosome and chromosome position
def extract_chr_and_position(vus_df: pd.DataFrame):
    # create new columns for chromosome and position and initialize it with empty strings
    vus_df.insert(1, 'Chr', "")
    vus_df.insert(2, 'Position', "")

    # make a copy of the dataframe to be able to iterate through it whilst modifying the original dataframe
    new_vus_df = vus_df.copy()

    # iterate through the dataframe and add entries to the chromosome and position column based on the locus column
    for index, row in new_vus_df.iterrows():
        # get the chromosome and the chromosome position from the locus
        locus = row['Locus']
        locus_arr = locus.split(':')

        vus_df.at[index, 'Chr'] = locus_arr[0].split('chr')[1]
        vus_df.at[index, 'Position'] = locus_arr[1]

    # remove locus column
    vus_df = vus_df.drop(columns=['Locus'])

    return vus_df


def preprocess_vus(vus_df: pd.DataFrame) -> InternalResponse:
    # exclude technical artifacts and CNVs
    vus_df = vus_df[vus_df['Classification'].str.contains('TECHNICAL_ARTIFACT', case=False, regex=True) == False]
    vus_df = vus_df[vus_df['Type'].str.contains('CNV|LONGDEL', case=False, regex=True) == False]

    # remove variant id column
    vus_df = vus_df.drop(columns=['Variant Id'])

    # remove flag type column
    vus_df = vus_df.drop(columns=['FlagType'])

    vus_df = extract_chr_and_position(vus_df)

    get_rsids_from_dbsnp_res = get_rsids_from_dbsnp(vus_df)

    if get_rsids_from_dbsnp_res.status != 200:
        current_app.logger.error(
            f'Get RSIDs from dbSNP query failed 500')
        return InternalResponse(None, 500)
    else:
        vus_df = get_rsids_from_dbsnp_res.data

        vus_df = handle_vus_with_multiple_genes(vus_df)

        retrieve_clinvar_variant_classifications_res = retrieve_clinvar_variant_classifications(vus_df)

        if retrieve_clinvar_variant_classifications_res.status != 200:
            current_app.logger.error(
                f'Retrieval of ClinVar variant classifications failed 500')
            return InternalResponse(None, 500)
        else:
            vus_df = retrieve_clinvar_variant_classifications_res.data

            return InternalResponse(vus_df, 200)


def handle_vus_file(file: FileStorage) -> Response:
    vus_df = pd.read_excel(file, header=2)  # TODO: check re header
    current_app.logger.info(f'Number of VUS found in file: {len(vus_df)}')

    # TODO: check - should I assume that inputted file only has VUS?

    preprocess_vus_res = preprocess_vus(vus_df)

    if preprocess_vus_res.status != 200:
        return Response({'isSuccess': False}, 200)
    else:
        vus_df = preprocess_vus_res.data

        # write dataframe to excel file
        # TODO: write to database
        vus_df.to_excel('final_vus.xlsx', index=False)

        return Response(json.dumps({'isSuccess': True}), 200)
