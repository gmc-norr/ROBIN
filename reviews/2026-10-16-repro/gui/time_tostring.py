import time, numpy as np, pandas as pd, hashlib, pickle, io
for n in (10_000, 50_000):
    rng = np.random.default_rng(0)
    df = pd.DataFrame({
        "read_id": [f"{i:08x}-aaaa-bbbb-cccc-{i:012d}" for i in range(n)],
        "col1": rng.choice(["chr1","chr2","chr7"], n), "col2": rng.integers(0, 2e8, n), "col3": rng.integers(0,2e8,n),
        "col4": rng.choice(["EGFR","SEPT14","FGFR3","TACC3"], n), "mapping_quality": rng.integers(0,60,n),
        "strand": rng.choice(["+","-"], n), "read_start": rng.integers(0,2e8,n), "read_end": rng.integers(0,2e8,n),
        "gene": rng.choice(["EGFR","SEPT14"], n), "chromosome": rng.choice(["chr1","chr7"], n), "start": rng.integers(0,2e8,n), "end": rng.integers(0,2e8,n),
    })
    b = pickle.dumps({"annotated_data": df})
    t=time.perf_counter(); pickle.loads(b); t1=time.perf_counter()-t
    t=time.perf_counter(); s=df.to_string(); hashlib.md5(s.encode()).hexdigest(); t2=time.perf_counter()-t
    print(f"rows={n}: unpickle {t1:.3f}s, to_string+md5 {t2:.2f}s, str MB {len(s)/1e6:.1f}")
