[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
[![DOI](https://img.shields.io/badge/doi-10.21468/SciPostPhysCodeb.64-darkblue.svg)](https://scipost.org/SciPostPhysCodeb.64)
[![arXiv](https://img.shields.io/badge/arXiv-2507.06082-b31b1b.svg)](https://arxiv.org/abs/2507.06082)
[![codecov](https://codecov.io/github/thomasbiekoetter/copa/graph/badge.svg?token=JC3O49KMQU)](https://codecov.io/github/thomasbiekoetter/copa)
[![last-commit](https://img.shields.io/github/last-commit/thomasbiekoetter/copa)](https://github.com/thomasbiekoetter/copa/commits/master)

# **copa – Chains for Optimization and Probabilistic Analysis**

**copa** is a lightweight and modular Fortran library for **Markov Chain Monte Carlo (MCMC)** sampling and probabilistic analysis.  
It provides parallel and serial ensemble samplers with convenient tools for storing and analyzing Markov chains. copa implements a parallel ensemble MCMC sampler based on the **Affine Invariant Ensemble Sampler** algorithm introduced by **Goodman & Weare (2010)** — *Communications in Applied Mathematics and Computational Science, 5(1), 65–80* ([DOI:10.2140/camcos.2010.5.65](https://doi.org/10.2140/camcos.2010.5.65)).

## **pycopa – Python Wrapper**

**pycopa** is a thin, user-friendly Python wrapper around the Fortran **copa** library. It exposes copa’s ensemble MCMC samplers. See the GitHub repository: [pycopa](https://github.com/thomasbiekoetter/pycopa)

---

## 📜 License and Citation

**copa** is licensed under the **GNU General Public License v3 (GPLv3)**.

If you use copa in academic work, please cite the accompanying paper on **evortran**:

> [arXiv:2507.06082]: Thomas Biekötter (IFT, Madrid), *evortran: a modern Fortran package for genetic algorithms with applications from LHC data fitting to LISA signal reconstruction*, [SciPost Phys. Codebases 64 (2026)]

```bibtex
@article{Biekotter:2025gkp,
    author = {Biek{\"o}tter, Thomas},
    title = "{evortran: A modern Fortran package for genetic algorithms with applications from LHC data fitting to LISA signal reconstruction}",
    eprint = "2507.06082",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "IFT-UAM/CSIC-25-76",
    doi = "10.21468/SciPostPhysCodeb.64",
    journal = "SciPost Phys. Codeb.",
    volume = "64",
    pages = "1",
    year = "2026"
}
```

## 🚀 Installation

You can build copa using the GNU gfortran or the Intel ifx compiler and the [Fortran Package Manager (fpm)](https://github.com/fortran-lang/fpm):

```bash
git clone https://gitlab.com/thomas.biekoetter/copa.git
fpm build --profile release
```
One can also build copa in debug mode, which contains additional compiler checks and runtime checks:
```bash
fpm build --profile debug
```

## 🧩 Example: Sampling the Rosenbrock Function

The example program `copa__test_rosenbrock` demonstrates how to use copa to sample the **2D Rosenbrock function**, a common benchmark for optimization and MCMC methods.

To run the example:

```bash
fpm test copa__test_rosenbrock
```

This will generate binary data files:

- `plots/rosenbrock/chains.npy` – the full chain data  
- `plots/rosenbrock/log_probs.npy` – the log-probability traces  

You can analyze these with **NumPy** and visualize results using tools like [corner.py](https://corner.readthedocs.io/) or [GetDist](https://getdist.readthedocs.io/en/latest/).

## ⚙️ Parallel vs Serial Sampler

copa includes both **parallel** and **serial** ensemble samplers:

- **Parallel (OpenMP)** – uses multiple threads to accelerate sampling:  
  ```fortran
  call run_parallel_sampler( &
      ndim, log_prior, log_like, &
      nsteps=nsteps, nthreads=nthreads, &
      ranges=ranges, &
      walkers=walkers, chains=chains, log_probs=log_probs)
  ```

- **Serial (single-threaded)** – same interface, simply omit the `nthreads` argument:  
  ```fortran
  call run_sampler( &
      ndim, log_prior, log_like, &
      nsteps=nsteps, &
      ranges=ranges, &
      walkers=walkers, chains=chains, log_probs=log_probs)
  ```

The number of threads can be controlled via OpenMP environment variables (e.g. `OMP_NUM_THREADS`) or at runtime within the code with the optional `nthreads` argument.

### 🧩 Input Arguments

Both `run_sampler` and `run_parallel_sampler` take a similar set of arguments.  
Some are **required**, while others are **optional**.

#### **Required Arguments**

- `ndim`  
  Integer — the number of free parameters (dimensions) in the model to be sampled.

- `log_prior`  
  User-defined subroutine with interface  
  ```fortran
  subroutine log_prior(theta, logp)
      real(wp), intent(in)  :: theta(:)
      real(wp), intent(out) :: logp
  end subroutine log_prior
  ```
  Defines the logarithm of the prior probability distribution.  
  Should return a large negative value (e.g. `-huge(1.0_wp)`) for invalid regions.

- `log_like`  
  User-defined subroutine with interface  
  ```fortran
  subroutine log_like(theta, logl)
      real(wp), intent(in)  :: theta(:)
      real(wp), intent(out) :: logl
  end subroutine log_like
  ```
  Defines the logarithm of the likelihood function.

#### **Optional Arguments**

- `nsteps`  
  Integer — the number of sampling steps to run for each walker (default: 1000).

- `nthreads`  
  Integer — the number of OpenMP threads to use (only for `run_parallel_sampler`).  
  If omitted, OpenMP will use the maximum number of available threads.

- `ranges`  
  Real array of shape `(2, ndim)` — defines lower and upper bounds for each parameter.  
  For example, to sample all parameters within zero and one:
  ```fortran
  ranges(1,:) = 0.0e0_wp   ! lower bounds
  ranges(2,:) = 1.0e0_wp   ! upper bounds
  ```

These arguments control the initialization, parallelization, and prior support region for the sampling run.

### Returned Arrays

Both samplers return their results through the arguments `walkers`, `chains`, and `log_probs`:

- `walkers([nthreads,] ndim, nwalkers)`  
  The current positions of the ensemble of walkers at the latest sampling step.  
  In the parallel sampler, the first dimension corresponds to the number of OpenMP threads.

- `chains([nthreads,] ndim, nwalkers, nsteps)`  
  The full sampling history of all walkers.  
  Each thread in the parallel version produces its own independent chain block along the fourth dimension.

- `log_probs([nthreads,] nwalkers, nsteps)`  
  The log-probability values corresponding to each sampled state.  
  Again, the first dimension is present only in the parallel sampler.

