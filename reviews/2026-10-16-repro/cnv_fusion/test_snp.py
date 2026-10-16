import sys, importlib.util
from pathlib import Path
spec = importlib.util.spec_from_file_location("snp_processing", "./src/robin/analysis/snp_processing.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
for v in ["Conflicting_interpretations_of_pathogenicity", "Conflicting_classifications_of_pathogenicity", "Benign", "Likely_pathogenic", "Pathogenic"]:
    print(f"{v:50s} -> is_pathogenic={m._is_pathogenic('CLNSIG='+v)}")
df = m.parse_vcf(Path(sys.argv[1]))
print(df[["CHROM","POS","Gene_Name" if "Gene_Name" in df.columns else "CHROM","is_pathogenic"]].to_string())
d = m.build_snp_display_data(Path(sys.argv[1]))
print("summary:", d["summary"])
