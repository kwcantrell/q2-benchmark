import biom
import matplotlib.pyplot as plt
import pandas as pd
import upsetplot


def _format_uplot(df: pd.DataFrame) -> pd.DataFrame:
    """Formats feature presence DataFrame for use in upsetplot

    Parameters
    ----------
    df : pd.DataFrame
        Pandas DataFrame of counts of each feature in each group

    Returns
    -------
    pd.DataFrame
        Pandas DataFrame formatted for use in upsetplot
    """
    df_binary = (df > 0.0).rename(columns=lambda x: x + ">0")
    df_up = pd.concat([df, df_binary], axis=1)
    df_up = df_up.set_index(list(df_binary.columns))
    return df_up


def _create_group_presence_df(
    table: biom.Table,
    metadata: pd.DataFrame,
    column: str
) -> pd.DataFrame:
    """Create DataFrame of each feature's presence in levels of a group

    Parameters
    ----------
    table : biom.Table
        Feature table
    metadata : pd.DataFrame
        Pandas DataFrame to use for group membership calculation
    column : str
        Column in metadata for group membership calculation

    Returns
    -------
    pd.DataFrame
        Pandas DataFrame of each feature and presence in each group
    """
    group_dfs = []
    for group, group_df in metadata.groupby(column):
        group_table = table.filter(group_df.index, inplace=False)
        non_zero_feats = group_table.sum("observation") > 0
        feats_to_keep = group_table.ids("observation")[non_zero_feats]
        group_table.filter(feats_to_keep, axis="observation")
        group_table_df = pd.DataFrame(
            group_table.sum(axis="observation"),
            group_table.ids("observation"),
            [group]
        )
        group_dfs.append(group_table_df)

    all_groups_df = pd.concat(group_dfs, axis=1).fillna(0)
    return all_groups_df


def plot(
    table: biom.Table,
    metadata: pd.DataFrame,
    column: str,
    save_loc: str
) -> None:
    """Create and save upset plot

    Parameters
    ----------
    table : biom.Table
        Feature table
    metadata : pd.DataFrame
        Pandas DataFrame to use for group membership calculation
    column : str
        Column in metadata for group membership calculation
    save_loc : str
        File to save plot
    """
    df = _create_group_presence_df(table, metadata, column)
    df = _format_uplot(df)

    fig = plt.figure()
    upset = upsetplot.UpSet(
        df,
        subset_size="count",
        element_size=100,
        show_counts=True
    )
    upset.plot(fig=fig)
    plt.savefig(save_loc)
