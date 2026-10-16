# Replicates the exact loops in src/robin/reporting/sections/fusion.py
# (_get_validated_fusion_pairs lines ~49-60, _get_gene_networks, _format_fusion_table ~249-262)
import pandas as pd
# clustered_data as returned by GUI _cluster_fusion_reads: sorted by supporting reads desc
clustered = pd.DataFrame([
  {"fusion_pair":"BRAF-KIAA1549","chr1":"chr7","chr2":"chr7","gene1_position":"140780000-140781000","gene2_position":"138860000-138861000","reads":25},
  {"fusion_pair":"BRAF-KIAA1549","chr1":"chr7","chr2":"chr7","gene1_position":"140790000-140791000","gene2_position":"138900000-138901000","reads":4},
  {"fusion_pair":"EGFR,EGFR-AS1-SEPTIN14","chr1":"chr7","chr2":"chr7","gene1_position":"55200000-55201000","gene2_position":"55800000-55801000","reads":9},
])
# summary path
validated_pairs=[]; seen=set()
for _, row in clustered.iterrows():
    s=row["fusion_pair"]
    if s and s not in seen:
        genes=[g.strip() for g in s.split("-") if g.strip()]
        if len(genes)>=2:
            key=tuple(sorted(genes))
            if key not in seen:
                validated_pairs.append(key); seen.add(key)
# connected components via union-find (equivalent to nx.connected_components on edges pair[0]-pair[1])
parent={}
def f(x):
    parent.setdefault(x,x)
    while parent[x]!=x: x=parent[x]
    return x
for p in validated_pairs: parent[f(p[0])]=f(p[1])
comps={}
for n in list(parent): comps.setdefault(f(n),[]).append(n)
print("validated_pairs:", validated_pairs)
print("summary networks:", [" - ".join(c) for c in comps.values()])
# detailed table path
details={}
for _, row in clustered.iterrows():
    s=row["fusion_pair"]
    genes=[g.strip() for g in s.split("-") if g.strip()]
    if len(genes)>=2:
        details["-".join(sorted(genes))]={"pos1":row["gene1_position"],"reads":int(row["reads"])}
print("PDF table rows:", details)
