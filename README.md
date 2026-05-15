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

`run_all.m` executes, in order: Pendulum (`Run_Pendulum_Batch` + `PendulumPlotting`), Tennis racket (`Run_Tennis_Full_Batch` + `Tennis_Plotting`), and Heavy Top (`Run_HeavyTop_batch` + `HeavyTop_plotting`). Each example is run in its own folder; allow up to 1 hour for a full run. You can still run any example individually using the README in that folder.

System used for development and testing:
- Operating system: Windows 10/11 (64‑bit)
- Machine type: Work laptop, multi‑core CPU
- MATLAB: R2024b (64‑bit), standard desktop installation
- Disk (`.mat` outputs only; not tracked in git): measured full batches are on the order of ~12 GB total. Recommended free space: **at least 4 GB** (Pendulum), **4 GB** (Tennis racket), **6 GB** (Heavy Top), or **13 GB** for all three together (includes modest headroom above typical runs).
