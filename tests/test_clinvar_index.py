"""VAR-3: a ClinVar refresh must not keep the tabix index of the old file.

A tabix index addresses records by 16 kb bin, so a stale index only *visibly*
hides the new data when the replacement record lands in a different bin from
the one it replaced. Position 100 -> 5,000,000 does; 100 -> 200 does not, and a
test built on the latter passes with or without the fix. Every VCF here is
hand-written and synthetic, and nothing touches the network.
"""

import os
import shutil
import time
from pathlib import Path

import pysam
import pytest

from robin.utils import clinvar_manager
from robin.utils.clinvar_manager import (
    CLINVAR_VCF_GZ_NAME,
    CLINVAR_VCF_GZ_TBI_NAME,
    _ensure_tabix_index,
)

_VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=249250621>\n"
    '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def _vcf_text(pos: int) -> str:
    return _VCF_HEADER + f"chr1\t{pos}\tsynthetic{pos}\tA\tG\t.\t.\tCLNSIG=Pathogenic\n"


def _bgzip_to(gz_path: Path, pos: int) -> None:
    """Write a one-record synthetic VCF to `gz_path`, bgzipped and unindexed."""
    gz_path.parent.mkdir(parents=True, exist_ok=True)
    plain_path = gz_path.parent / f"{gz_path.name}.source.vcf"
    plain_path.write_text(_vcf_text(pos))
    pysam.tabix_compress(str(plain_path), str(gz_path), force=True)
    plain_path.unlink()


def _write_indexed_vcf(directory: Path, *, pos: int) -> tuple[Path, Path]:
    """A bgzipped ClinVar-shaped VCF plus a matching, correct tabix index."""
    gz_path = directory / CLINVAR_VCF_GZ_NAME
    _bgzip_to(gz_path, pos)
    pysam.tabix_index(str(gz_path), preset="vcf", force=True, keep_original=True)
    return gz_path, directory / CLINVAR_VCF_GZ_TBI_NAME


def _overwrite_gz_only(gz_path: Path, *, pos: int) -> None:
    """Swap a newer .gz in under the surviving .tbi, as a ClinVar refresh does.

    The .tbi is also backdated an hour. That is the real situation -- the index
    was built for the file that has just been replaced -- and it is what makes
    "the index was rewritten" assertable by mtime at all: this machine stores
    mtime to ~1 ms, so a rebuild landing in the same millisecond as the original
    leaves `st_mtime_ns` untouched -- measured here at 122 of 300 runs without
    the backdating, and 0 of 300 with it.
    """
    staged_gz_path = gz_path.parent / "refresh.vcf.gz"
    _bgzip_to(staged_gz_path, pos)
    shutil.copyfile(staged_gz_path, gz_path)
    staged_gz_path.unlink()

    tbi_path = gz_path.parent / CLINVAR_VCF_GZ_TBI_NAME
    if tbi_path.exists():
        an_hour_ago = int((time.time() - 3600) * 1e9)
        os.utime(tbi_path, ns=(an_hour_ago, an_hour_ago))


def _refuse_download(url, target_path, **kwargs):
    raise AssertionError(f"no download was expected, but one was attempted: {url}")


def test_a_replaced_vcf_gets_a_rebuilt_index(tmp_path):
    """VAR-3: the old .tbi survives the swap, so the new records are invisible."""
    gz, tbi = _write_indexed_vcf(tmp_path, pos=100)
    _overwrite_gz_only(gz, pos=5_000_000)  # the .tbi is now stale
    _ensure_tabix_index(gz, tbi, force=True)
    with pysam.TabixFile(str(gz)) as tf:
        assert [r.split("\t")[1] for r in tf.fetch("chr1", 4_999_000, 5_001_000)] == ["5000000"]


def test_force_actually_rewrites_the_index_file(tmp_path):
    """Guards against a `force` parameter that is accepted and then ignored."""
    gz, tbi = _write_indexed_vcf(tmp_path, pos=100)
    _overwrite_gz_only(gz, pos=5_000_000)
    before = tbi.stat().st_mtime_ns
    _ensure_tabix_index(gz, tbi, force=True)
    assert tbi.stat().st_mtime_ns != before


def test_without_force_a_good_index_is_not_rebuilt(tmp_path):
    gz, tbi = _write_indexed_vcf(tmp_path, pos=100)
    before = tbi.stat().st_mtime_ns
    _ensure_tabix_index(gz, tbi)
    assert tbi.stat().st_mtime_ns == before


def test_the_stale_index_really_hides_the_new_record(tmp_path):
    """Pins the bug itself: this is what users saw -- zero results, exit code 0."""
    gz, tbi = _write_indexed_vcf(tmp_path, pos=100)
    _overwrite_gz_only(gz, pos=5_000_000)
    with pysam.TabixFile(str(gz)) as tf:
        assert list(tf.fetch("chr1", 4_999_000, 5_001_000)) == []


def test_ensure_clinvar_files_reindexes_a_freshly_downloaded_gz(tmp_path, monkeypatch):
    """The `ensure_clinvar_files` call site: a new .gz arrives beside a stale .tbi."""
    resources_dir = tmp_path / "resources"
    gz, tbi = _write_indexed_vcf(resources_dir, pos=100)
    gz.unlink()  # the .gz is gone; the index built for it survives

    def fake_download(url, target_path, *, timeout_s=600):
        _bgzip_to(Path(target_path), 5_000_000)

    monkeypatch.setattr(clinvar_manager, "_download_url_to_file", fake_download)

    clinvar_manager.ensure_clinvar_files(
        resources_dir=resources_dir, url="http://clinvar.invalid/clinvar.vcf.gz"
    )

    with pysam.TabixFile(str(gz)) as tf:
        assert [r.split("\t")[1] for r in tf.fetch("chr1", 4_999_000, 5_001_000)] == ["5000000"]


def test_update_clinvar_if_newer_reindexes_after_the_swap(tmp_path, monkeypatch):
    """The `update_clinvar_if_newer` call site, immediately after the replace()."""
    resources_dir = tmp_path / "resources"
    gz, tbi = _write_indexed_vcf(resources_dir, pos=100)

    def fake_download(url, target_path, *, timeout_s=600):
        _bgzip_to(Path(target_path), 5_000_000)

    monkeypatch.setattr(clinvar_manager, "_download_url_to_file", fake_download)
    monkeypatch.setattr(
        clinvar_manager,
        "_get_remote_last_modified",
        lambda url, **kwargs: time.time() + 86_400,
    )

    assert (
        clinvar_manager.update_clinvar_if_newer(
            resources_dir=resources_dir, url="http://clinvar.invalid/clinvar.vcf.gz"
        )
        is True
    )

    with pysam.TabixFile(str(gz)) as tf:
        assert [r.split("\t")[1] for r in tf.fetch("chr1", 4_999_000, 5_001_000)] == ["5000000"]


def test_an_already_present_gz_is_not_reindexed_on_start_up(tmp_path, monkeypatch):
    """Guard -- passes before and after the fix.

    It pins the half of the `ensure_clinvar_files` change that has no
    observable behaviour of its own: when the .gz was *not* just downloaded,
    `force` must stay False so a good index is not rebuilt on every start-up.
    That cannot be shown through the index, because the only way to reach this
    call without a download is with the .tbi already missing, where `force`
    makes no difference -- hence the recorded argument.
    """
    resources_dir = tmp_path / "resources"
    gz, tbi = _write_indexed_vcf(resources_dir, pos=100)
    tbi.unlink()  # the .gz is fine; only the index is missing

    real_ensure_tabix_index = clinvar_manager._ensure_tabix_index
    recorded_force_values: list[object] = []

    def recording_ensure(gz_path, tbi_path, **kwargs):
        recorded_force_values.append(kwargs.get("force", False))
        return real_ensure_tabix_index(gz_path, tbi_path, **kwargs)

    monkeypatch.setattr(clinvar_manager, "_download_url_to_file", _refuse_download)
    monkeypatch.setattr(clinvar_manager, "_ensure_tabix_index", recording_ensure)

    clinvar_manager.ensure_clinvar_files(
        resources_dir=resources_dir, url="http://clinvar.invalid/clinvar.vcf.gz"
    )

    assert recorded_force_values == [False]
    assert tbi.exists()


def _raise_instead_of_indexing(*args, **kwargs):
    raise OSError("simulated tabix indexing failure (disk full)")


def test_a_failed_forced_rebuild_leaves_no_index_behind(tmp_path, monkeypatch):
    """I-1: a forced rebuild that fails must take the stale .tbi with it.

    The rebuild raising is correct and loud. What is not correct is leaving the
    index of the *previous* release on disk: it passes `_file_ok`, so the next
    start-up trusts it. Nothing here asserts on mtimes -- this machine stores
    mtime to ~1 ms and a rebuild often lands in the same millisecond -- only on
    whether the file is there.
    """
    gz, tbi = _write_indexed_vcf(tmp_path, pos=100)
    _overwrite_gz_only(gz, pos=5_000_000)  # the .tbi now indexes the old file
    assert tbi.exists()

    monkeypatch.setattr(pysam, "tabix_index", _raise_instead_of_indexing)
    monkeypatch.setattr(clinvar_manager.shutil, "which", lambda name: None)

    with pytest.raises(RuntimeError):
        _ensure_tabix_index(gz, tbi, force=True)

    assert not tbi.exists()


def test_a_failed_rebuild_is_not_trusted_on_the_next_start_up(tmp_path, monkeypatch):
    """I-1, the property that matters: VAR-3 must not come back one run later.

    An indexing failure is usually transient (a full disk, a truncated
    download). The run after it must notice the index is missing and rebuild,
    not find a stale .tbi passing `_file_ok` and return early -- which would
    leave every region query returning no records and exit code 0, with the
    only error already in a previous run's log.
    """
    resources_dir = tmp_path / "resources"
    gz, tbi = _write_indexed_vcf(resources_dir, pos=100)
    _overwrite_gz_only(gz, pos=5_000_000)

    real_tabix_index = pysam.tabix_index
    real_which = clinvar_manager.shutil.which
    monkeypatch.setattr(pysam, "tabix_index", _raise_instead_of_indexing)
    monkeypatch.setattr(clinvar_manager.shutil, "which", lambda name: None)
    with pytest.raises(RuntimeError):
        _ensure_tabix_index(gz, tbi, force=True)

    # Next start-up, indexing working again and no download to be had.
    monkeypatch.setattr(pysam, "tabix_index", real_tabix_index)
    monkeypatch.setattr(clinvar_manager.shutil, "which", real_which)

    real_ensure_tabix_index = clinvar_manager._ensure_tabix_index
    recorded_force_values: list[object] = []

    def recording_ensure(gz_path, tbi_path, **kwargs):
        recorded_force_values.append(kwargs.get("force", False))
        return real_ensure_tabix_index(gz_path, tbi_path, **kwargs)

    monkeypatch.setattr(clinvar_manager, "_download_url_to_file", _refuse_download)
    monkeypatch.setattr(clinvar_manager, "_ensure_tabix_index", recording_ensure)

    clinvar_manager.ensure_clinvar_files(
        resources_dir=resources_dir, url="http://clinvar.invalid/clinvar.vcf.gz"
    )

    # Reached at all -- not short-circuited by the early "both files OK" return.
    assert recorded_force_values == [False]
    with pysam.TabixFile(str(gz)) as tf:
        assert [r.split("\t")[1] for r in tf.fetch("chr1", 4_999_000, 5_001_000)] == ["5000000"]
