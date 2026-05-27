# LGIMsAccuracy
This study conducts a systematic comparison of six methods: conventional second-order differential equations and Newton-Euler angular velocity-based equations, each implemented with both Euler parameters and Cartesian rotation vectors, alongside corresponding Lie group integration (LGIM) variants using RKMK integration.

Three problems are being analyzed (each one is in separate folder):
1. **Pendulum** - constrained pendulum exhibiting smooth, large rotations
<br> Read *READMEPendulum.txt* to make comparison on Pendulum problem

2. **Tennis racket phenomena** also known as the intermediate axis theorem or Dzhanibekov effect,
describes the unstable rotational behavior of a rigid body about its intermediate principal axis.
<br> Read *READMETennisRacket.txt* to make comparison on Tennis racket phenomena

3. **Heavy Top** - classical benchmark problem in multibody system dynamics that combines constrained motion with fast rotations and rapid rotational changes.
<br> Read *READMEHeavyTop.txt* to make comparison on Heavy Top problem

## Run all examples at once

From the repository root, you can run every batch and plotting script in sequence:

```matlab
cd('<path-to>/LGIMsAccuracy')   % repository root
run_all
```

`run_all.m` executes, in order: Pendulum (`Run_Pendulum_Batch` + `PendulumPlotting`), Tennis racket (`Run_Tennis_Full_Batch` + `Tennis_Plotting`), and Heavy Top (`Run_HeavyTop_batch` + `HeavyTop_plotting`). Each example is run in its own folder. You can still run any example individually using the README in that folder.

## Computational time (`executionTime`)

Each integrator records **`executionTime`** in its `.mat` file. Plotting scripts read this field for “computation time vs time step” curves.

**What is timed:** only the **`custom_rk4` integration** (the main time-stepping loop). Initial setup, constraint/energy post-processing, and saving `.mat` files are **not** included. This keeps cross-method timing comparable.

**Repeats:** before integration, each method script checks workspace variable **`n_timing_repeats`** (set by the batch runner, or default `1` if you run a single script by hand). The integration runs that many times; **`executionTime`** is the **median** of those samples (robust to occasional slow runs from background load). Trajectories are taken from the **last** repeat only—repeats affect timing statistics, not accuracy results.

**Batch default:** `n_timing_repeats = 1` in `Run_Pendulum_Batch.m`, `Run_Tennis_Full_Batch.m`, and `Run_HeavyTop_batch.m`. For publication-style timing plots with less scatter, set e.g. **`n_timing_repeats = 5`** in the relevant batch file. In practice, increasing repeats scales integration cost by about **×N** (saving and post-processing are not repeated), so a full batch with 5 repeats can take on the order of **several hours** on a laptop.

**Typical batch runtimes** (Windows laptop, MATLAB R2024b, `n_timing_repeats = 1`; plotting extra):

| Example        | Batch (approx.) | Notes                          |
|----------------|-----------------|--------------------------------|
| Pendulum       | ~10 min         | See `READMEPendulum.txt`       |
| Tennis racket  | ~10 min         | See `READMETennisRacket.txt`   |
| Heavy Top      | 30–45 min       | See `READMEHeavyTop.txt`       |
| All three      | ~1 h            | `run_all` + plotting per folder |

System used for development and testing:
- Operating system: Windows 10/11 (64‑bit)
- Machine type: Work laptop, multi‑core CPU
- MATLAB: R2024b (64‑bit), standard desktop installation
- Disk (`.mat` outputs only; not tracked in git): measured full batches are on the order of ~12 GB total. Recommended free space: **at least 4 GB** (Pendulum), **4 GB** (Tennis racket), **6 GB** (Heavy Top), or **13 GB** for all three together (includes modest headroom above typical runs).
