import numpy as np
import pandas as pd
import corner
import matplotlib.pyplot as plt
from getdist import MCSamples, plots
import glob
import re


ndim = 2
nwalkers = 50
nsteps = 10000
nthreads = 4

burn_in = nsteps / 10

files = glob.glob("chains_*.npy")
files.sort(key=lambda f: [int(x) for x in re.findall(r'\d+', f)])
samples = [np.fromfile(f, dtype=np.float64).reshape((nsteps, ndim)) for f in files]

files = glob.glob("log_probs_*.npy")
files.sort(key=lambda f: [int(x) for x in re.findall(r'\d+', f)])
loglikes = [np.fromfile(f, dtype=np.float64).reshape((nsteps, )) for f in files]

labels = []
for i in range(0, ndim):
    l = r'\theta_' + str(i + 1)
    labels.append(l)

samples = MCSamples(
    samples=samples,
    labels=labels,
    ignore_rows=burn_in,
    loglikes=loglikes,
    names=['theta1', 'theta2'])

g = plots.get_subplot_plotter()

g.triangle_plot(
    [samples],
    shaded=True,
    filled=True,
    title_limit=1)

g.export("bla.pdf")

r_minus_1 = samples.getGelmanRubin()
print("Gelman-Rubin R-1 = ", r_minus_1)

marge = samples.getMargeStats()
print(marge.parWithName('param1').mean)
