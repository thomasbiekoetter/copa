import numpy as np
import matplotlib.pyplot as plt


ndim = 2
nwalkers = 50
nsteps = 10000
nthreads = 8

size = 100000

samples = np.fromfile('chains.npy', dtype=np.float64)
samples = samples.reshape((nwalkers * nsteps * nthreads, ndim))
samples = samples[-size:-1, :]

log_probs = np.fromfile('log_probs.npy', dtype=np.float64)
log_probs = log_probs[-size:-1]

fig, ax = plt.subplots()

sc = ax.scatter(
    samples[:, 0],
    samples[:, 1],
    c=log_probs,
    s=4,
    rasterized=True)
fig.colorbar(sc)

plt.savefig('log_probs.pdf')


