import ast, time, logging, hashlib, numpy as np, pandas as pd
src = open("./src/robin/gui/components/fusion.py").read()
tree = ast.parse(src)
want = {"_extract_fusion_breakpoints","_cluster_breakpoints","_cluster_positions","_position_in_cluster",
        "_positions_overlap","_positions_close","_validate_fusion_breakpoints","_create_data_hash"}
mod = ast.Module(body=[n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in want], type_ignores=[])
ns = {"pd": pd, "np": np, "logging": logging, "hashlib": hashlib, "List": list, "Dict": dict, "Any": object, "Optional": object}
exec(compile(mod, "fusion_funcs", "exec"), ns)
rng = np.random.default_rng(1)
for nreads in (100, 200, 400):
    rows = []
    pairs = [("EGFR","SEPT14"),("FGFR3","TACC3"),("KIAA1549","BRAF"),("NTRK2","X1")]
    for i in range(nreads):
        g1, g2 = pairs[i % len(pairs)]
        base1, base2 = rng.integers(1e6, 1e8), rng.integers(1e6, 1e8)
        if i % 10 < 7: base1, base2 = 5_000_000 + int(rng.integers(0,50)), 9_000_000 + int(rng.integers(0,50))  # recurrent breakpoint
        for g, b in ((g1, base1), (g2, base2)):
            rows.append({"read_id": f"r{i}", "col4": g, "reference_id": "chr7", "reference_start": int(b),
                         "reference_end": int(b) + 500, "strand": "+", "mapping_quality": 60, "mapping_span": 500})
    df = pd.DataFrame(rows)
    t = time.perf_counter(); v = ns["_validate_fusion_breakpoints"](df, 4, 100); t1 = time.perf_counter() - t
    t = time.perf_counter(); ns["_create_data_hash"]({"annotated_data": df, "candidate_count": 1}); t2 = time.perf_counter() - t
    print(f"reads={nreads} rows={len(df)}: validate {t1:.2f}s, hash {t2:.2f}s (x2 files, every 30s per client, on event loop)")
