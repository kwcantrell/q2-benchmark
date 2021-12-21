from pkg_resources import resource_filename

import biom
import pandas as pd
import pytest
from qiime2 import Artifact

from q2benchmark import upset

DATA_DIR = resource_filename("q2benchmark", "tests/test_data")


@pytest.fixture
def load_data():
    tbl = Artifact.load(f"{DATA_DIR}/extraction_test_round_3_biom_lod.qza")
    tbl = tbl.view(biom.Table)

    md = pd.read_table(
        f"{DATA_DIR}/metadata_12201_round3_qiitaIDs_2020.08.12_qiime2.txt",
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

    tax = pd.read_table(
        f"{DATA_DIR}/taxonomy.tsv",
        sep="\t",
        index_col=0
    ).loc[tbl.ids(axis="observation")]
    return tbl, md, tax


def test_create_group_presence_df(load_data):
    tbl, md, _ = load_data
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
    tbl, md, _ = load_data
    df = upset._create_group_presence_df(
        table=tbl,
        metadata=md,
        column="extraction_kit_bead_beating"
    )
    u_df = upset._format_uplot(df)

    exp_cols = set(md["extraction_kit_bead_beating"])
    act_cols = set(u_df.columns)
    assert exp_cols == act_cols


def test_taxsplit(load_data):
    tbl, md, tax = load_data
    df = upset._create_group_presence_df(
        table=tbl,
        metadata=md,
        column="extraction_kit_bead_beating"
    )
    exp_rows = [3, 49, 141, 238, 319, 874, 628]  # k, p, c, o, f, g, s
    for level, tax_exp_rows in enumerate(exp_rows, start=1):
        df = upset._taxsplit(df, tax, level)
        assert df.shape == (tax_exp_rows, 3)


def test_taxsplit_level_high(load_data):
    tbl, md, tax = load_data
    df = upset._create_group_presence_df(
        table=tbl,
        metadata=md,
        column="extraction_kit_bead_beating"
    )
    with pytest.raises(ValueError) as excinfo:
        upset._taxsplit(df, tax, 8)
    exp_error_msg = "Level must be less than total number of levels!"
    assert str(excinfo.value) == exp_error_msg


def test_taxsplit_level_low(load_data):
    tbl, md, tax = load_data
    df = upset._create_group_presence_df(
        table=tbl,
        metadata=md,
        column="extraction_kit_bead_beating"
    )

    for i in [-1, 0]:
        with pytest.raises(ValueError) as excinfo:
            upset._taxsplit(df, tax, 0)
        exp_error_msg = "Level cannot be <= 0!"
        assert str(excinfo.value) == exp_error_msg


def test_plot_no_tax_file(tmpdir, load_data):
    tbl, md, tax = load_data
    out_file = tmpdir.mkdir("sub").join("test.png")

    with pytest.raises(ValueError) as excinfo:
        upset.plot(tbl, md, "extraction_kit_bead_beating", out_file,
                   tax_level=1)
    exp_error_msg = "Must provide taxonomy!"
    assert str(excinfo.value) == exp_error_msg


def test_plot(tmpdir, load_data):
    tbl, md, _ = load_data
    out_file = tmpdir.mkdir("sub").join("test.png")
    upset.plot(tbl, md, "extraction_kit_bead_beating", out_file)
