import numpy as np
rng=np.random.default_rng(0)
# chr7 at 1 Mb bins (~160 bins), ploidy ~2 with noise; EGFR focal amplification (2 bins at ~30 copies)
v=rng.normal(2,0.3,160); v[55:57]=30
mean,std=v.mean(),v.std()
print("per-chromosome plot ylim = (0, %.2f); amplified bins at 30 are %s" % (mean+2*std, "CLIPPED (not drawn)" if 30>mean+2*std else "visible"))
# genome-wide: threshold = mean + 5*std over all bins
g=np.concatenate([rng.normal(2,0.3,3000), v])
print("genome-wide ylim upper = %.2f" % (g.mean()+5*g.std()))
