from pkg_resources import resource_filename

import biom
import pandas as pd
import pytest
from qiime2 import Artifact

from q2benchmark import upset


@pytest.fixture
def load_data():
    data_dir = resource_filename("q2benchmark", "tests/test_data")
    tbl = Artifact.load(f"{data_dir}/extraction_test_round_3_biom_lod.qza")
    tbl = tbl.view(biom.Table)
    md = pd.read_table(
        f"{data_dir}/metadata_12201_round3_qiitaIDs_2020.08.12_qiime2.txt",
        sep="\t",
        index_col=0,
        skiprows=[1],
    )
    md["extraction_kit_bead_beating"] = (
        md.extraction_kit + "\n" + md.bead_beating
    )
    samps_to_keep = set(tbl.ids()).intersection(md.index)

    tbl.filter(samps_to_keep)
    md = md.loc[samps_to_keep]
    return tbl, md


def test_create_group_presence_df(load_data):
    tbl, md = load_data
    df = upset._create_group_presence_df(
        table=tbl,
        metadata=md,
        column="extraction_kit_bead_beating"
    )
    assert df.shape == (37933, 3)

    # Order of features not guaranteed to be the same
    exp_rows = set(df.index)
    act_rows = set(tbl.ids(axis="observation"))
    assert exp_rows == act_rows

    exp_cols = set(md["extraction_kit_bead_beating"])
    act_cols = set(df.columns)
    assert exp_cols == act_cols


def test_format_uplot(load_data):
    tbl, md = load_data
    df = upset._create_group_presence_df(
        table=tbl,
        metadata=md,
        column="extraction_kit_bead_beating"
    )
    u_df = upset._format_uplot(df)

    exp_cols = set(md["extraction_kit_bead_beating"])
    act_cols = set(u_df.columns)
    assert exp_cols == act_cols


def test_plot(tmpdir, load_data):
    tbl, md = load_data
    out_file = tmpdir.mkdir("sub").join("test.png")
    upset.plot(tbl, md, "extraction_kit_bead_beating", out_file)
