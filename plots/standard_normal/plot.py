import numpy as np
import pandas as pd
import corner
import matplotlib.pyplot as plt


ndim = 3
nwalkers = 50
nsteps = 10000

burn_in = 100

samples = np.fromfile('chains.npy', dtype=np.float64)
samples = samples.reshape((nwalkers*nsteps, ndim))
samples = samples[burn_in:-1,:]

figure = corner.corner(
    samples,
    bins=40,
    smooth=1,
    labels=["theta1", "theta2", "theta3"],  # replace with your param names
    show_titles=True,
    title_fmt=".2f",
    levels=(0.68, 0.95),  # confidence levels
    range=[(-3, 3), (-3, 3), (-3, 3)],
    plot_datapoints=False)

plt.savefig('corner.pdf')
