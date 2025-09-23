import numpy as np
import pandas as pd
import corner
import matplotlib.pyplot as plt


ndim = 3
nwalkers = 20
nsteps = 100000

burn_in = 10000

samples = np.fromfile('chains.npy', dtype=np.float64)
samples = samples.reshape((nwalkers*nsteps, ndim))
samples = samples[burn_in:-1,:]

figure = corner.corner(
    samples,
    labels=["theta1", "theta2", "theta3"],  # replace with your param names
    show_titles=True,
    title_fmt=".2f",
#   quantiles=[0.16, 0.5, 0.84],  # 1-sigma intervals
    levels=(0.68, 0.95),  # confidence levels
    plot_datapoints=False)

db = pd.read_csv("evortran_best.csv").columns.astype(float).to_numpy()
axs = figure.get_axes()
ax = axs[3]
ax.scatter(
    db[0], db[1],
    s=400,
    marker='*',
    zorder=1000,
    color='magenta')
ax = axs[6]
ax.scatter(
    db[0], db[2],
    s=400,
    marker='*',
    zorder=1000,
    color='magenta')
ax = axs[7]
ax.scatter(
    db[1], db[2],
    s=200,
    marker='*',
    zorder=1000,
    color='magenta')


plt.savefig('corner.pdf')
