import numpy as np
import pandas as pd
import corner
import matplotlib.pyplot as plt
from getdist import MCSamples, plots



ndim = 2
nwalkers = 50
nsteps = 10000
nthreads = 4

burn_in = nsteps / 10

samples_1 = np.fromfile("chains.npy_1_1.npy", dtype=np.float64)
samples_2 = np.fromfile("chains.npy_2_1.npy", dtype=np.float64)
samples_3 = np.fromfile("chains.npy_3_1.npy", dtype=np.float64)
samples_4 = np.fromfile("chains.npy_4_1.npy", dtype=np.float64)
samples_5 = np.fromfile("chains.npy_5_1.npy", dtype=np.float64)
samples_6 = np.fromfile("chains.npy_6_1.npy", dtype=np.float64)
samples_7 = np.fromfile("chains.npy_7_1.npy", dtype=np.float64)
samples_8 = np.fromfile("chains.npy_8_1.npy", dtype=np.float64)

log_probs = np.fromfile('log_probs.npy', dtype=np.float64)
log_probs = log_probs[-size:-1]

samples = [
    samples_1,
    samples_2,
    samples_3,
    samples_4,
    samples_5,
    samples_6,
    samples_7,
    samples_8]

for i, s in enumerate(samples):
    s = s.reshape((nsteps, ndim))
    samples[i] = s

labels = []
for i in range(0, ndim):
    l = r'$\theta_' + str(i + 1) + r'$'
    labels.append(l)

samples = MCSamples(
    samples=samples,
    labels=labels,
    ignore_rows=burn_in)
#   loglikes=log_probs,
#   names=['T', 'K', 'bH'],
#   labels=['T', 'K', 'b/H'])

g = plots.get_subplot_plotter()

g.triangle_plot(
    [samples],
    filled=True)

g.export("bla.pdf")

quit()



# samples = np.fromfile('chains.npy', dtype=np.float64)
# samples = samples.reshape((nwalkers * nsteps * nthreads, ndim))
# samples = samples[burn_in:-1,:]

ranges = []
dlim = 0
ulim = 5
for i in range(0, ndim):
    ranges.append((dlim, ulim))



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
