import glob, re
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


ndim = 2
nwalkers = 100
nsteps = 1000
nensembles = 4

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

c = np.concatenate(chains, axis=0)
l = -np.concatenate(loglikes)

fig, ax = plt.subplots()
sc = ax.scatter(
    c[:, 0],
    c[:, 1],
    c=l,
    norm=LogNorm(vmin=l.min(), vmax=l.max()),
    s=4,
    rasterized=True)
fig.colorbar(sc)
plt.savefig('log_probs.pdf')


