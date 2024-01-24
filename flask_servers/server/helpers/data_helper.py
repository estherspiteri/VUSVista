import pandas as pd
from typing import List


def convert_df_to_list(df: pd.DataFrame) -> List:
    # convert dataframe to list
    final_list = []
    for index, row in df.iterrows():
        final_list.append(row.to_dict())
    print(final_list)
    return final_list


def prep_vus_df_for_react(vus_df: pd.DataFrame) -> pd.DataFrame:
    # to match React camelCase syntax
    new_vus_df = pd.DataFrame()
    new_vus_df['chromosome'] = vus_df['Chr']
    new_vus_df['chromosomePosition'] = vus_df['Position']
    new_vus_df['gene'] = vus_df['Gene']
    new_vus_df['type'] = vus_df['Type']
    new_vus_df['genotype'] = vus_df['Genotype']
    new_vus_df['refAllele'] = vus_df['Reference']
    new_vus_df['observedAllele'] = vus_df['Observed Allele']
    new_vus_df['classification'] = vus_df['Classification']
    new_vus_df['rsid'] = vus_df['RSID']
    new_vus_df['rsidDbsnpVerified'] = vus_df['RSID dbSNP verified']
    new_vus_df['rsidDbsnpErrorMsgs'] = vus_df['RSID dbSNP errorMsgs']
    new_vus_df['clinvarErrorMsg'] = vus_df['Clinvar error msg']
    new_vus_df['clinvarClassification'] = vus_df['Clinvar classification']
    new_vus_df['clinvarClassificationLastEval'] = vus_df['Clinvar classification last eval']
    new_vus_df['clinvarClassificationReviewStatus'] = vus_df['Clinvar classification review status']
    new_vus_df['clinvarCanonicalSpdi'] = vus_df['Clinvar canonical spdi']
    new_vus_df['clinvarUid'] = vus_df['Clinvar uid']

    return new_vus_df
