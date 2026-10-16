import sys, tempfile, pathlib, pickle
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent / "stubs"))
sys.path.insert(0, "./src")
import numpy as np
from robin.gui.components import summary as s
d = pathlib.Path(tempfile.mkdtemp())
# Simulate CNV outputs after 37 analysis iterations
(d / "cnv_analysis_counter.txt").write_text("37")
np.save(d / "CNV_dict.npy", {"bin_width": 1000000, "variance": 0.12})
pickle.dump("Male", open(d / "XYestimate.pkl", "wb"))
print("CNV summary:", s._extract_cnv_data(d))
# Coverage thresholds
import pandas as pd
for depth in (12.0, 25.0):
    pd.DataFrame({"chrom":["chr1"],"startpos":[1],"endpos":[1000],"name":["G"],"bases":[depth*1000]}).to_csv(d/"bed_coverage_main.csv", index=False)
    print(depth, "summary quality:", s._extract_coverage_data(d)["quality"])
# MGMT with final_mgmt.csv
pd.DataFrame({"average":[12.3],"pred":[45.6],"status":["Methylated"]}).to_csv(d/"final_mgmt.csv", index=False)
(d/"final_mgmt.bed").write_text("chr10\t1\t2\n")
print("MGMT summary:", s._extract_mgmt_data(d))
