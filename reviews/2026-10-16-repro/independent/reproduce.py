"""Reproduce additional review findings without importing ROBIN's GUI stack.

Run from any directory with Python 3.12 and numpy, pandas, pysam and pyarrow:
    python reviews/2026-10-16-repro/independent/reproduce.py -v

These are bug demonstrations: passing assertions confirm the current incorrect
behavior. They should be inverted when converted into regression tests for fixes.
Definitions are compiled unchanged from the working tree using ast. Real BAMs,
indexes, coverage commands, parquet files and file locks are used where relevant.
External services/R/classifiers and unrelated orchestration are stubbed explicitly.
All generated files, including the harmless shell marker, stay in temporary dirs.
"""

from __future__ import annotations

import array
import ast
import contextlib
import dataclasses
import fcntl
import gc
import glob
import io
import json
import logging
import os
from pathlib import Path
import pickle
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import types
import unittest
from unittest.mock import patch
import uuid

import numpy as np
import pandas as pd
import pysam
from dateutil import parser


ROOT = Path(__file__).resolve().parents[3]
LOG = logging.getLogger("independent-review")


def definitions(relative_path, names=None, **injected):
    """Load source definitions and literal constants, skipping package imports."""
    path = ROOT / relative_path
    tree = ast.parse(path.read_text(), filename=str(path))
    module = types.ModuleType(f"review_{path.stem}_{uuid.uuid4().hex}")
    module.__file__ = str(path)
    sys.modules[module.__name__] = module
    env = module.__dict__
    env.update(
        os=os, sys=sys, time=time, re=re, json=json, logging=logging,
        tempfile=tempfile, gc=gc, glob=glob, shutil=shutil, uuid=uuid,
        fcntl=fcntl, pd=pd, np=np, pysam=pysam, Path=Path, StringIO=io.StringIO,
        dataclass=dataclasses.dataclass, field=dataclasses.field, parser=parser,
        pickle=pickle, logger=LOG,
        resources=types.SimpleNamespace(__file__=str(ROOT / "src/robin/resources/__init__.py")),
        get_job_logger=lambda *args: LOG,
    )
    selected = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if names is None or node.name in names:
                selected.append(node)
        elif isinstance(node, (ast.Assign, ast.AnnAssign)):
            try:
                ast.literal_eval(node.value)
            except (ValueError, TypeError):
                continue
            selected.append(node)
    future = ast.ImportFrom(module="__future__", names=[ast.alias(name="annotations")], level=0)
    code = ast.fix_missing_locations(ast.Module(body=[future, *selected], type_ignores=[]))
    exec(compile(code, str(path), "exec", dont_inherit=True), env)
    env.update(injected)
    return module


def make_bam(path, records, *, length=1000, sample="sample", index=True):
    """records: (name, reference_start, sequence, flag, optional tags)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": length}, {"SN": "chr2", "LN": length}],
        "RG": [{"ID": "run", "LB": sample, "DS": "runid=run basecall_model=test"}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for name, pos, seq, flag, *rest in records:
            read = pysam.AlignedSegment(bam.header)
            read.query_name = name
            read.query_sequence = seq
            read.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            read.flag = flag
            read.reference_id = 0
            read.reference_start = pos
            read.mapping_quality = 60
            read.cigarstring = f"{len(seq)}M"
            read.set_tag("RG", "run")
            for key, value in (rest[0] if rest else {}).items():
                read.set_tag(key, value)
            bam.write(read)
    if index:
        pysam.index(str(path))
    return path


def target_module():
    return definitions("src/robin/analysis/target_analysis.py", {
        "FileLock", "TargetMetadata", "_load_bed_regions", "run_bedtools",
        "get_covdfs", "get_read_counts_per_target", "run_bedmerge", "TargetAnalysis",
    })


def make_target(module, directory, bed):
    worker = module.TargetAnalysis(work_dir=str(directory), target_panel="rCNS2", batch_size=1)
    worker.bedfile = str(bed)
    # These fixtures have no dynamic master BED; keep unrelated imports out.
    worker._get_master_bed_path = lambda sample_id: None
    return worker


class AdditionalFindings(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="robin-independent-review-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.output = io.StringIO()
        self.stdout = contextlib.redirect_stdout(self.output)
        self.stdout.__enter__()
        self.addCleanup(self.stdout.__exit__, None, None, None)

    def test_ir01_reference_only_sample_is_exported_as_pathogenic(self):
        module = definitions("src/robin/analysis/lightweight_gene_analysis.py", {
            "GeneReport", "LightweightGeneAnalysis", "is_supported", "needs_review", "resolve_sample_folder",
        })
        vcf = self.root / "clinvar.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t11\t123\tC\tT\t.\t.\tCLNSIG=Pathogenic;GENEINFO=TEST:1;CLNDN=test\n"
        )
        bam = make_bam(self.root / "normal.bam", [
            (f"r{i}", 0, "C" * 50, 0) for i in range(20)
        ])
        bed = self.root / "targets.bed"
        bed.write_text("chr1\t0\t50\tTEST\n")
        worker = module.LightweightGeneAnalysis(str(self.root), str(vcf))
        reports = worker.analyze_bed_regions(str(bam), str(bed))
        output = worker.export_json(reports, str(self.root / "result.json"))
        data = json.loads(Path(output).read_text())
        variant = data["genes"]["TEST"]["variants"][0]
        self.assertEqual(variant["evidence_analysis"]["alternate_support"], 0)
        self.assertEqual(variant["evidence_analysis"]["evidence_level"], "no_support")
        self.assertEqual(data["genes"]["TEST"]["summary"]["pathogenic_status"], "present_with_good_coverage")
        self.assertEqual(data["metadata"]["genes_with_pathogenic_variants"], 1)

    def test_ir02_missing_sample_bam_falls_back_to_another_sample(self):
        module = definitions("src/robin/analysis/lightweight_gene_analysis.py", {
            "GeneReport", "LightweightGeneAnalysis", "is_supported", "needs_review", "resolve_sample_folder", "lightweight_gene_analysis_handler",
        })
        wrong_bam = make_bam(self.root / "B/igv/igv_ready.bam", [
            (f"B-read{i}", 0, "T" * 30, 0) for i in range(20)
        ])
        sample_a = self.root / "A"
        sample_a.mkdir()
        make_bam(sample_a / "target.bam", [(f"A-read{i}", 0, "C" * 30, 0) for i in range(20)])
        bed = sample_a / "targets_exceeding_threshold.bed"
        bed.write_text("chr1\t0\t30\tTEST\n")
        vcf = self.root / "clinvar.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t11\t123\tC\tT\t.\t.\tCLNSIG=Pathogenic;GENEINFO=TEST:1;CLNDN=test\n"
        )
        observed = {}

        class RecordInputs(module.LightweightGeneAnalysis):
            def analyze_bed_regions(self, bam, targets):
                observed.update(bam=bam, bed=targets)
                return super().analyze_bed_regions(bam, targets)

            def export_json(self, reports, output):
                observed["output"] = output
                return super().export_json(reports, output)

        module.LightweightGeneAnalysis = RecordInputs
        context = types.SimpleNamespace(
            filepath=str(sample_a), metadata={"bed_path": str(bed), "clinvar_path": str(vcf)},
            get_sample_id=lambda: "A", add_result=lambda *args: None,
            add_error=lambda *args: self.fail(str(args)),
        )
        module.lightweight_gene_analysis_handler(
            types.SimpleNamespace(job_id=1, job_type="lightweight_gene_analysis", context=context),
            work_dir=str(self.root),
        )
        self.assertEqual(observed["bam"], str(wrong_bam))
        self.assertEqual(observed["bed"], str(bed))
        self.assertEqual(observed["output"], str(sample_a / "lightweight_gene_analysis_results.json"))
        data = json.loads(Path(observed["output"]).read_text())
        evidence = data["genes"]["TEST"]["variants"][0]["evidence_analysis"]
        self.assertEqual(evidence["alternate_support"], 20)  # from B, written under A

    def test_ir03_global_depth_depends_on_bam_partition(self):
        module = target_module()
        bed = self.root / "panel.bed"
        bed.write_text("chr1\t0\t100\tTEST\n")
        records = [(f"r{i}", 0, "C" * 100, 0) for i in range(10)]
        single = make_bam(self.root / "all.bam", records, length=100)
        halves = [make_bam(self.root / f"half{i}.bam", rows, length=100)
                  for i, rows in enumerate((records[:5], records[5:]))]
        results = []
        for label, bams in (("single", [single]), ("split", halves)):
            worker = make_target(module, self.root / label, bed)
            for bam in bams:
                metadata, _ = worker.process_file_with_staging(str(bam), {"sample_id": "sample"})
                self.assertIsNone(metadata.error_message)
                result = worker.accumulate_staged_files("sample", force=True)
                self.assertEqual(result["status"], "success")
            # Two equal-sized contigs; chr2 has zero reads. True global depth = 5x.
            results.append(result["coverage"])
        self.assertEqual(results, [0.5, 1.0])
        cov, _ = module.get_covdfs(str(single), str(bed))
        self.assertEqual(cov.loc[cov["#rname"] == "chr1", "meandepth"].iloc[0], 10)

    def test_ir04_concurrent_target_accumulations_overwrite_a_batch(self):
        module = target_module()
        bed = self.root / "panel.bed"
        bed.write_text("chr1\t0\t100\tTEST\n")
        work = self.root / "work"
        first = make_target(module, work, bed)
        second = make_target(module, work, bed)
        # Suppress BAM creation to isolate the coverage read/modify/write race.
        first.callthreshold = second.callthreshold = float("inf")
        bams = [make_bam(self.root / f"input{i}.bam", [(f"r{i}", 0, "C" * 100, 0)], length=100)
                for i in range(2)]
        barrier = threading.Barrier(2)
        real_counts = module.get_read_counts_per_target

        def wait_after_existing_state_was_read(*args):
            barrier.wait(timeout=10)
            return real_counts(*args)

        module.get_read_counts_per_target = wait_after_existing_state_was_read
        results = []
        metadata, _ = first.process_file_with_staging(str(bams[0]), {"sample_id": "sample"})
        self.assertIsNone(metadata.error_message)
        thread = threading.Thread(target=lambda: results.append(first.accumulate_staged_files("sample", force=True)))
        thread.start()
        # Wait until first accumulation has claimed its files and read old state.
        deadline = time.monotonic() + 10
        while barrier.n_waiting != 1:
            if time.monotonic() > deadline:
                self.fail("First accumulation did not reach the synchronization point")
            time.sleep(0.01)
        metadata, _ = second.process_file_with_staging(str(bams[1]), {"sample_id": "sample"})
        self.assertIsNone(metadata.error_message)
        results.append(second.accumulate_staged_files("sample", force=True))
        thread.join(timeout=10)
        self.assertFalse(thread.is_alive())
        self.assertEqual([r["status"] for r in results], ["success", "success"])
        coverage = pd.read_csv(work / "sample/bed_coverage_main.csv")
        self.assertEqual(int(coverage["bases"].sum()), 100)  # expected 200
        self.assertEqual(first._get_pending_count("sample"), 0)

    def test_ir05_implicit_canonical_mm_calls_are_omitted(self):
        module = definitions("src/robin/analysis/utilities/matkit.py", {
            "map_read_to_ref_positions", "process_bam_counts_improved",
        })
        # Each read covers the same two CpGs. MM '.' declares unlisted C bases
        # canonical; one read explicitly modifies CpG1, the other CpG2.
        bam = make_bam(self.root / "implicit.bam", [
            ("r1", 0, "CGCG", 0, {"MM": "C+m.,0;", "ML": array.array("B", [255])}),
            ("r2", 0, "CGCG", 0, {"MM": "C+m.,1;", "ML": array.array("B", [255])}),
        ])
        counts, _, _ = module.process_bam_counts_improved(str(bam), combine_mods=True)
        for pos in (0, 2):
            row = counts[("chr1", pos, "+", "m")]
            self.assertEqual(row[module.COUNT_IDX_MOD], 1)
            self.assertEqual(row[module.COUNT_IDX_CANONICAL], 0)  # expected 1
            self.assertEqual(row[module.COUNT_IDX_TOTAL], 1)  # expected 2

    def test_ir06_low_quality_supplementary_call_depends_on_record_order(self):
        module = definitions("src/robin/analysis/utilities/matkit.py", {
            "map_read_to_ref_positions", "process_bam_counts_improved",
        })
        tags = {"MM": "C+m,0;", "ML": array.array("B", [255]), "qs": 5.0}
        # Both files are coordinate sorted; the primary switches genomic location.
        before = make_bam(self.root / "supp_first.bam", [
            ("r", 0, "CG", 2048, tags), ("r", 10, "CG", 0, tags),
        ])
        after = make_bam(self.root / "primary_first.bam", [
            ("r", 0, "CG", 0, tags), ("r", 10, "CG", 2048, tags),
        ])
        before_counts, _, _ = module.process_bam_counts_improved(str(before), combine_mods=True)
        after_counts, _, _ = module.process_bam_counts_improved(str(after), combine_mods=True)
        self.assertEqual(len(before_counts), 1)
        self.assertEqual(len(after_counts), 0)

    def test_ir07_supplementary_alignment_is_counted_as_another_read(self):
        module = definitions("src/robin/analysis/bam_preprocessor.py", {
            "get_rg_tags_from_bam", "process_bam_reads",
        }, _BARCODE_PATTERN=re.compile(r"_barcode(\d{1,2})$"), _LJ_BAM_THREADS=0)
        bam = make_bam(self.root / "pass.bam", [
            ("one-read", 0, "CGCG", 0), ("one-read", 10, "CGCG", 2048),
        ])
        data = module.process_bam_reads(str(bam))
        self.assertEqual(data["mapped_reads"], 2)  # one physical read
        self.assertEqual(data["yield_tracking"], 8)  # four sequenced bases

    def test_ir08_mgmt_sample_path_is_executed_by_the_shell(self):
        module = definitions("src/robin/analysis/mgmt_analysis.py", {"run_matkit", "convert_to_mixed_delim"})
        sample = self.root / "sample$(touch SHELL_MARKER)"
        sample.mkdir()
        hv = self.root / "hv"
        (hv / "bin").mkdir(parents=True)
        (hv / "bin/mgmt_pred_v0.3.R").touch()

        def sort_stub(**kwargs):
            Path(kwargs["output_bam"]).touch()
            return True

        def matkit_stub(input_bam, output):
            Path(output).write_text("\t".join(["0"] * 18) + "\n")

        module.safely_sort_and_index_bam = sort_stub
        # Run the actual shell with a harmless Rscript substitute, in a temp cwd.
        module.subprocess = types.SimpleNamespace(run=lambda cmd, **kwargs: subprocess.run(
            cmd, **kwargs, cwd=self.root, env={**os.environ, "PATH": str(self.root)}))
        (self.root / "Rscript").write_text("#!/bin/sh\nexit 0\n")
        (self.root / "Rscript").chmod(0o755)
        (self.root / "touch").symlink_to(shutil.which("touch"))
        stub = types.ModuleType("robin.analysis.utilities.matkit")
        stub.run_matkit = matkit_stub
        with patch.dict(sys.modules, {stub.__name__: stub}):
            self.assertTrue(module.run_matkit(str(self.root), "unused.bam", str(hv), str(sample), 1))
        self.assertTrue((self.root / "SHELL_MARKER").exists())

    def test_ir09_repeated_region_name_only_pileups_first_interval(self):
        module = definitions("src/robin/analysis/lightweight_gene_analysis.py", {
            "GeneReport", "LightweightGeneAnalysis", "is_supported", "needs_review", "resolve_sample_folder",
        })
        vcf = self.root / "clinvar.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t11\t1\tC\tT\t.\t.\tCLNSIG=Pathogenic;GENEINFO=TEST:1\n"
            "1\t111\t2\tC\tT\t.\t.\tCLNSIG=Pathogenic;GENEINFO=TEST:1\n"
        )
        bed = self.root / "two_exons.bed"
        bed.write_text("chr1\t0\t30\tTEST\nchr1\t100\t130\tTEST\n")
        bam = make_bam(self.root / "two_exons.bam", [
            (f"r{i}", pos, "T" * 30, 0)
            for pos in (0, 100) for i in range(5)
        ])
        worker = module.LightweightGeneAnalysis(str(self.root), str(vcf))
        reports = worker.analyze_bed_regions(str(bam), str(bed))
        self.assertEqual([v["position"] for v in reports["TEST"].pathogenic_snps_details], [11])

    def test_ir10_cnv_master_bed_excludes_the_selected_panel_without_fusion(self):
        master = definitions("src/robin/analysis/master_bed_generator.py")
        cnv = definitions("src/robin/analysis/cnv_analysis.py", {
            "find_significant_regions", "generate_bed_files", "save_cnv_files",
        })
        sample_dir = self.root / "sample"
        sample_dir.mkdir()
        (sample_dir / "master.csv").write_text("analysis_panel\ndemo\n")
        panel = self.root / "panel.bed"
        panel.write_text("chr1\t0\t100\tTEST\n")
        master._get_target_bed_path = lambda name: str(panel)
        fusion_stub = types.ModuleType("robin.analysis.fusion_work")
        fusion_stub._load_fusion_metadata = lambda *args: None
        fusion_stub._load_analysis_counter = lambda *args: 1
        master_stub = types.ModuleType("robin.analysis.master_bed_generator")
        master_stub._try_get_target_panel_from_fusion_metadata = master._try_get_target_panel_from_fusion_metadata
        calls = []

        def generate_synchronously(**kwargs):
            calls.append(kwargs)
            return master.generate_master_bed(**kwargs)

        # Execute background generation synchronously for a deterministic fixture.
        master_stub.generate_master_bed_async = generate_synchronously
        with patch.dict(sys.modules, {
            fusion_stub.__name__: fusion_stub, master_stub.__name__: master_stub,
        }):
            cnv.save_cnv_files(
                str(sample_dir), 1, {"chr2": [2.0] * 10}, {"chr2": [2.0] * 10},
                {"chr2": [0.0] * 10}, 100, 0.0,
                [{"chromosome": "chr2", "start": 200, "end": 400}],
                "Unknown", {}, LOG, reference=str(self.root / "missing-reference.fa"),
                generate_master_bed=True,
            )
            self.assertIsNone(calls[0]["target_panel"])
            target = target_module().TargetAnalysis(work_dir=str(self.root), target_panel="rCNS2")
            selected = target._get_master_bed_path("sample")
            self.assertIsNotNone(selected)
            self.assertNotIn("chr1", Path(selected).read_text())
            # The existence cache also prevents a subsequent correct request
            # from repairing the same counter's incomplete master BED.
            master.generate_master_bed("sample", str(self.root), 1, target_panel="demo")
            self.assertNotIn("chr1", Path(selected).read_text())
            bam = make_bam(self.root / "panel-read.bam", [("r", 0, "C" * 100, 0)])
            filtered = self.root / "filtered.bam"
            target_module().run_bedtools(str(bam), selected, str(filtered))
            with pysam.AlignmentFile(str(filtered), "rb") as result:
                self.assertEqual(result.count(until_eof=True), 0)

    def test_ir11_prethreshold_bams_never_enter_the_variant_calling_input(self):
        module = target_module()
        bed = self.root / "panel.bed"
        bed.write_text("chr1\t0\t1000\tTEST\n")
        worker = make_target(module, self.root / "work", bed)
        for i in range(3):
            bam = make_bam(self.root / f"early{i}.bam", [(f"r{i}", 0, "C" * 50, 0)])
            metadata, _ = worker.process_file_with_staging(str(bam), {"sample_id": "sample"})
            self.assertIsNone(metadata.error_message)
            self.assertEqual(worker.accumulate_staged_files("sample", force=True)["status"], "success")
        batches = list((self.root / "work/sample").glob("batch_*.bam"))
        names = []
        for path in batches:
            with pysam.AlignmentFile(str(path), "rb") as bam:
                names.extend(read.query_name for read in bam.fetch(until_eof=True))
        self.assertEqual(names, ["r2"])  # r0 and r1 are never revisited
        self.assertEqual(worker._get_pending_count("sample"), 0)


if __name__ == "__main__":
    logging.basicConfig(level=logging.CRITICAL)
    unittest.main()
