import numpy as np
import pandas as pd
import corner
import matplotlib.pyplot as plt


ndim = 2
nwalkers = 50
nsteps = 10000
nthreads = 8

burn_in = 1000

samples = np.fromfile('chains.npy', dtype=np.float64)
samples = samples.reshape((nwalkers * nsteps * nthreads, ndim))
samples = samples[burn_in:-1,:]

ranges = []
dlim = -1
ulim = 3.5
for i in range(0, ndim):
    ranges.append((dlim, ulim))

labels = []
for i in range(0, ndim):
    l = r'$\theta_' + str(i + 1) + r'$'
    labels.append(l)

figure = corner.corner(
    samples,
    bins=40,
    smooth=1,
    labels=labels,
    show_titles=True,
    title_fmt=".2f",
    levels=(0.68, 0.95),
    range=ranges,
#   flat=True,
    plot_datapoints=False,
)

plt.savefig('corner.pdf')
