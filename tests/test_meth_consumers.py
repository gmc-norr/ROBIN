"""Locked rule 6 (plan a95935d): each classifier's input carries one row per
panel CpG, with the cumulative parquet's percent on that classifier's scale.

The parquet is built by the real merge from per-BAM rows. CpG 1300's strands
disagree (1 read at 100% on +, 6 reads at 0% on -): combined it is 14.3%
methylated. Before strands were combined, NanoDX saw only the + row (100%),
Sturgeon saw two opposite votes (no call), and random forest and MNP-Flex got
two rows per CpG. Classifiers, models and R scripts aren't run: each input is
caught where it would be handed over. Synthetic rows only.
"""

import importlib
import sys
import types

import numpy as np
import pandas as pd
import pytest

from meth_fixtures import write_per_bam_parquet
from robin.analysis.temp_utilities import merge_modkit_files

# pos -> (+ strand n_mod, n_canonical), (- strand n_mod, n_canonical)
SITES = {900: ((3, 1), (2, 0)), 1300: ((1, 0), (0, 6)), 1700: ((2, 2), (1, 1)), 2100: ((6, 0), None)}
PERCENT = {900: 100 * 5 / 6, 1300: 100 * 1 / 7, 1700: 50.0, 2100: 100.0}
COVERAGE = {900: 6, 1300: 7, 1700: 6, 2100: 6}


@pytest.fixture
def parquet(tmp_path):
    panel = tmp_path / "parquet_filter.txt"
    panel.write_text("chr start end IlmnID\n" + "".join(f"chr1\t{p}\t{p + 2}\tcg{p}\n" for p in SITES))
    rows = []
    for pos, (plus, minus) in SITES.items():
        rows.append((pos, "+", *plus))
        if minus:
            rows.append((pos + 1, "-", *minus))
    per_bam = write_per_bam_parquet(tmp_path / "b.parquet", rows)
    (tmp_path / "work" / "S").mkdir(parents=True)
    out = str(tmp_path / "work" / "S" / "S.parquet")
    merge_modkit_files([per_bam], out, out, str(panel), "S", str(tmp_path / "work"), None, 1)
    return out


@pytest.fixture
def mb(monkeypatch):
    """merge_bedmethyl, which imports sturgeon.utils at the top. Where sturgeon
    isn't installed, a stand-in module is used for this test only."""
    try:
        import sturgeon.utils  # noqa: F401
    except ImportError:
        stub, utils = types.ModuleType("sturgeon"), types.ModuleType("sturgeon.utils")

        def _not_used(*args, **kwargs):
            raise NotImplementedError("sturgeon stand-in")

        utils.read_probes_file = _not_used
        stub.utils = utils
        monkeypatch.setitem(sys.modules, "sturgeon", stub)
        monkeypatch.setitem(sys.modules, "sturgeon.utils", utils)
    name = "robin.analysis.utilities.merge_bedmethyl"
    had = name in sys.modules
    yield importlib.import_module(name)
    if not had:
        sys.modules.pop(name, None)


def test_nanodx_gets_one_call_per_cpg_from_both_strands(parquet, mb, monkeypatch, tmp_path):
    from robin.analysis import nanodx_analysis as nd

    captured = {}

    def classify(modelfile, test_df):
        captured["df"] = test_df.copy()
        return np.array([1.0]), np.array(["class"]), len(test_df)

    monkeypatch.setattr(nd, "load_modkit_data", mb.load_modkit_data)
    monkeypatch.setattr(nd, "collapse_minimal_bedmethyl", mb.collapse_minimal_bedmethyl)
    monkeypatch.setattr(nd, "_direct_classification", classify)
    analysis = nd.NanodxAnalysis.__new__(nd.NanodxAnalysis)
    analysis.model, analysis.modelfile = "Capper_et_al_NN_v2.pkl", "unused"
    analysis.work_dir, analysis.storefile = str(tmp_path), "NanoDX_scores.csv"
    analysis.cpgs = pd.DataFrame({0: ["chr1"] * 4, 1: list(SITES), 2: [p + 2 for p in SITES], 3: [f"cg{p}" for p in SITES]})

    analysis.process_parquet_file(parquet, "S")
    df = captured["df"]
    assert sorted(df["start_pos"]) == sorted(SITES)
    calls = dict(zip(df["start_pos"], df["methylation_call"]))
    assert calls == {p: (1 if PERCENT[p] >= 60 else -1) for p in SITES}


def test_sturgeon_gets_one_vote_per_cpg(parquet, mb, tmp_path):
    probes = tmp_path / "probes.bed"
    probes.write_text("chr start end probe_name\n" + "".join(f"chr1 {p} {p + 1} cg{p}\n" for p in SITES))
    result = mb.modkit_pileup_file_to_bed(mb.load_modkit_data(parquet), str(tmp_path / "pileup"), str(probes))
    # Sturgeon's rule per vote: above 70% methylated, below 30% unmethylated.
    expected = {f"cg{p}": int(PERCENT[p] > 70) for p in SITES if PERCENT[p] > 70 or PERCENT[p] < 30}
    assert dict(zip(result["probe_id"], result["methylation_call"])) == expected


def test_random_forest_gets_one_row_per_cpg_at_its_percent(parquet, mb):
    collapsed = mb.collapse_bedmethyl(mb.reconstruct_full_bedmethyl_data(mb.load_modkit_data(parquet)))
    assert sorted(collapsed["start_pos"]) == sorted(SITES)
    for pos, score in zip(collapsed["start_pos"], collapsed["score"]):
        assert score == pytest.approx(PERCENT[pos], abs=1e-4)
        assert score <= 100


def test_mnp_flex_gets_one_row_per_cpg_with_its_coverage(parquet, tmp_path):
    from robin.analysis.utilities.matkit import reconstruct_full_bedmethyl_for_mnpflex
    from robin.analysis.utilities.mnp_flex import APIClient

    bed = tmp_path / "S.mnpflex.bed"
    reconstruct_full_bedmethyl_for_mnpflex(parquet).to_csv(bed, sep="\t", index=False, header=False)
    windows = tmp_path / "windows.bed"
    windows.write_text("chr start end coverage methylation_percentage IlmnID\n"
                       + "".join(f"chr1 {p} {p + 2} 0 0 cg{p}\n" for p in SITES))
    out = tmp_path / "subset.bed"
    APIClient(base_url="https://example.invalid", verify_ssl=False).process_streaming(str(windows), str(bed), str(out))
    rows = [line.split() for line in out.read_text().splitlines()[1:]]
    assert sorted(int(r[1]) for r in rows) == sorted(SITES)
    for r in rows:
        pos = int(r[1])
        assert int(r[3]) == COVERAGE[pos]
        assert float(r[4]) == pytest.approx(PERCENT[pos], abs=0.005)
