import numpy as np
import matplotlib.pyplot as plt
import glob
import re


ndim = 2
nwalkers = 50
nsteps = 10000
nthreads = 4

files = glob.glob("chains_*.npy")
files.sort(key=lambda f: [int(x) for x in re.findall(r'\d+', f)])
samples = [np.fromfile(f, dtype=np.float64).reshape((nsteps, ndim)) for f in files]
samples = np.concatenate(samples, axis=0)

files = glob.glob("log_probs_*.npy")
files.sort(key=lambda f: [int(x) for x in re.findall(r'\d+', f)])
loglikes = [np.fromfile(f, dtype=np.float64).reshape((nsteps, )) for f in files]
loglikes = np.concatenate(loglikes, axis=0)

fig, ax = plt.subplots()

sc = ax.scatter(
    samples[:, 0],
    samples[:, 1],
    c=loglikes,
    s=4,
    rasterized=True)
fig.colorbar(sc)

plt.savefig('log_probs.pdf')


