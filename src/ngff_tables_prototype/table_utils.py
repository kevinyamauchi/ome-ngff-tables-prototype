from typing import Any, Dict, List, Optional

from anndata import AnnData
import pandas as pd


def df_to_anndata(
    df: pd.DataFrame,
    dense_columns: List[str],
    var: Optional[pd.DataFrame] = None
) -> AnnData:
    """Create the anndata object from the example csv
    Parameters
    ----------
    df : pd.DataFrame
        The table to convert to an AnnData object
    dense_columns : List[str]
        The list of column names to remove from the dataframe and
        store in the dense matrix (X)
    var : Optional[pd.DataFrame]
        The annotations for the columns of the dense matrix X that
        get stored in AnnData.var. Must be ordered to match
        dense_columns.
    Returns
    -------
    ann_obj : AnnData
        The converted AnnData object
    """
    # get the dense array
    dense_array = df[dense_columns].to_numpy()

    # drop the dense array from the table
    obs = df.drop(dense_columns, axis='columns')

    # create the AnnData object
    ann_data_kwargs = {'X': dense_array, 'obs': obs}
    if var is not None:
        ann_data_kwargs['var'] = var
    ann_obj = AnnData(**ann_data_kwargs)

    return ann_obj
