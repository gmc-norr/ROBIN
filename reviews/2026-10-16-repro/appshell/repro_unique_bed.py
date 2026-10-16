"""Extract cli._generate_unique_gene_bed / _validate_bed_file verbatim and run on sample panels."""
import ast
import sys
import tempfile
from pathlib import Path
import types

SRC = Path("./src/robin/cli.py").read_text()
tree = ast.parse(SRC)
wanted = {"_generate_unique_gene_bed", "_validate_bed_file"}
mod = ast.Module(body=[n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in wanted], type_ignores=[])
click_stub = types.SimpleNamespace(echo=lambda *a, **k: print("click.echo:", *a))
ns = {"click": click_stub, "Path": Path, "Tuple": tuple, "List": list}
from typing import Tuple, List
ns.update(Tuple=Tuple, List=List)
exec(compile(mod, "cli_extract", "exec"), ns)

tmp = Path(tempfile.mkdtemp(dir=sys.argv[1]))

# Case 1: typical multi-exon panel, 4 columns
bed4 = tmp / "panel4.bed"
bed4.write_text(
    "chr7\t55019017\t55019365\tEGFR\n"
    "chr7\t55142285\t55142437\tEGFR\n"
    "chr7\t55181292\t55181478\tEGFR\n"
    "chr17\t7669608\t7669690\tTP53\n"
    "chr17\t7673700\t7673837\tTP53\n"
)
print("validate 4col:", ns["_validate_bed_file"](bed4))
out4 = tmp / "out4.bed"
ns["_generate_unique_gene_bed"](bed4, out4)
print("--- output (4col input, 5 regions in):")
print(out4.read_text())

# Case 2: 8-column BED (passes validation since >=6 cols with valid score/strand)
bed8 = tmp / "panel8.bed"
bed8.write_text(
    "chr7\t55019017\t55019365\tEGFR\t0\t+\t55019017\t55019365\n"
    "chr17\t7669608\t7669690\tTP53\t0\t-\t7669608\t7669690\n"
)
print("validate 8col:", ns["_validate_bed_file"](bed8))
out8 = tmp / "out8.bed"
ok = ns["_generate_unique_gene_bed"](bed8, out8)
print("--- output (8col input):", ok)
print(out8.read_text() if out8.exists() else "<none>")
