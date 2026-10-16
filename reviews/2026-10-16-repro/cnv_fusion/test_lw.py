import sys, types, importlib.util, logging
# stub heavy deps not installed in the scratch venv
sys.modules['pysam'] = types.ModuleType('pysam')
lc = types.ModuleType('robin.logging_config'); lc.get_job_logger = lambda *a, **k: logging.getLogger('x')
sys.modules['robin'] = types.ModuleType('robin'); sys.modules['robin.logging_config'] = lc
spec = importlib.util.spec_from_file_location("lw", "./src/robin/analysis/lightweight_gene_analysis.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
obj = m.LightweightGeneAnalysis.__new__(m.LightweightGeneAnalysis)
obj.clinvar_variants = {}; obj.gene_variants = {}; obj.logger = logging.getLogger('x')
with open(sys.argv[1]) as f:
    obj._parse_clinvar_vcf(f)
print("clinvar_variants stored:", {k: (v['reference'], v['alternate'], v['clinical_significance']) for k, v in obj.clinvar_variants.items()})
print("gene_variants IDH1 count:", len(obj.gene_variants.get('IDH1', [])))
