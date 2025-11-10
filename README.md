# **copa – Chains for Optimization and Probabilistic Analysis**

**copa** is a lightweight and modular Fortran library for **Markov Chain Monte Carlo (MCMC)** sampling and probabilistic analysis.  
It provides parallel and serial ensemble samplers with convenient tools for storing and analyzing Markov chains. copa implements a parallel ensemble MCMC sampler based on the **Affine Invariant Ensemble Sampler** algorithm introduced by **Goodman & Weare (2010)** — *Communications in Applied Mathematics and Computational Science, 5(1), 65–80* ([DOI:10.2140/camcos.2010.5.65](https://doi.org/10.2140/camcos.2010.5.65)).


---

## 🚀 Installation

You can build copa using the gfortran compiler and the [Fortran Package Manager (fpm)](https://github.com/fortran-lang/fpm):

```bash
git clone https://gitlab.com/thomas.biekoetter/copa.git
cd copa
source exports_run_gfortran.sh
fpm build
```
Alternatively, one can build copa using the intel ifx compiler:
```bash
source exports_run_ifx.sh
fpm build
```

## 🧩 Example: Sampling the Rosenbrock Function

The example program `copa__test_rosenbrock` demonstrates how to use copa to sample the **2D Rosenbrock function**, a common benchmark for optimization and MCMC methods.

To run the example:

```bash
fpm test copa__test_rosenbrock
```

This will generate binary output files:

- `plots/rosenbrock/chains.npy` – the full chain data  
- `plots/rosenbrock/log_probs.npy` – the log-probability traces  

You can analyze these with **NumPy** and visualize results using tools like [corner.py](https://corner.readthedocs.io/).

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

- `walkers(ndim, nwalkers [, nthreads])`  
  The current positions of the ensemble of walkers at the latest sampling step.  
  In the parallel sampler, the third dimension corresponds to the number of OpenMP threads.

- `chains(ndim, nwalkers, nsteps [, nthreads])`  
  The full sampling history of all walkers.  
  Each thread in the parallel version produces its own independent chain block along the fourth dimension.

- `log_probs(nwalkers, nsteps [, nthreads])`  
  The log-probability values corresponding to each sampled state.  
  Again, the optional fourth dimension is present only in the parallel sampler.

## 📜 License and Citation

**copa** is licensed under the **GNU General Public License v3 (GPLv3)**.

If you use copa in academic work, please cite the accompanying paper on **evortran**:

> T. Biekoetter, *Evortran: Evolutionary Optimization and Random Sampling in Modern Fortran*, [arXiv:2507.06082](https://arxiv.org/abs/2507.06082)



