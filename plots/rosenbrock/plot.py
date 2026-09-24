import glob, re
import numpy as np
from collections import defaultdict
from getdist import MCSamples, plots

nsteps = 1000
ndim = 2
burn_in = nsteps // 5

names  = ['p1', 'p2']
labels = [r'p_1', r'p_2']

def idx(f): # 'chains_3_7.npy' -> (3, 7)
    iw, ie = map(int, re.findall(r'\d+', f))
    return iw, ie

chain_files = {idx(f): f for f in glob.glob("chains_*_*.npy")}
lp_files    = {idx(f): f for f in glob.glob("log_probs_*_*.npy")}

by_ensemble = defaultdict(list)
for (iw, ie) in chain_files:
    by_ensemble[ie].append(iw)

ensemble_samples = {}
for ie, walkers in sorted(by_ensemble.items()):
    chains   = [np.fromfile(chain_files[(iw, ie)], dtype=np.float64).reshape(nsteps, ndim)
                for iw in walkers]
    loglikes = [np.fromfile(lp_files[(iw, ie)], dtype=np.float64).reshape(nsteps)
                for iw in walkers]
    ensemble_samples[ie] = MCSamples(
        samples=chains, # list -> walkers pooled, kept separate for convergence
        loglikes=[-lp for lp in loglikes], # = -log(posterior) in getdist
        names=names,
        labels=labels,
        label=f'Ensemble {ie}',
        ignore_rows=burn_in) # applied per walker

samples = [ensemble_samples[ie] for ie in sorted(ensemble_samples)]

g = plots.get_subplot_plotter()

g.triangle_plot(
    samples,
    filled=True,
    title_limit=1)

g.export("corner.pdf")
