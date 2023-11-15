import pandas as pd
from typing import List


def convert_df_to_list(df: pd.DataFrame) -> List:
    # convert dataframe to list
    final_list = []
    for index, row in df.iterrows():
        final_list.append(row.to_dict())

    return final_list
