READMEPendulum

The pendulum example is a constrained multibody benchmark with smooth, large-amplitude rotation.
It is used here to compare six integration formulations (conventional and Newton–Euler style,
with Euler parameters and Cartesian rotation vectors, plus LGIM variants) against a
high-accuracy numerical reference.

System used for development and testing:
- Operating system: Windows 10/11 (64-bit)
- Machine type: Work laptop, multi-core CPU
- MATLAB: R2024b (64-bit), standard desktop installation
- Disk: Full batch `.mat` outputs are typically ~3 GB; allow **at least 4 GB** free (headroom).

Instructions

1) Copy all files from this `Pendulum` directory into a single folder on your local machine
   (or clone the GitHub repository and work directly in the `Pendulum` folder).
   Alternatively, from the repository root, run `run_all` to execute Pendulum, Tennis racket,
   and Heavy Top in one go (see the main `README.md`).
2) Open MATLAB R2024b (or a compatible recent version) and set the current folder to this
   `Pendulum` directory.
3) Run `Run_Pendulum_Batch` (the function defined in `Run_Pendulum_Batch.m`).
   This script will:
   - First run `ThetaRefode45.m` to build a high-accuracy reference and save it as
     `ExactRef_HighAccuracy.mat`.
   - Then simulate six formulations over one shared time-step list (milliseconds):
     `0.2, 0.5, 1, 2, 5, 10, 20, 50, 100` ms for Classic EP, Classic Cartesian, NE EP,
     wrapped NE Cartesian, Cartesian LGIM, and EP LGIM.
   - For each run, write a `.mat` file with the naming patterns expected by `PendulumPlotting.m`
     (`EP_dt_*.mat`, `Cart_dt_*.mat`, `EP_NE_dt_*.mat`, `wrCart_NE_dt_*.mat`,
     `wrLGIM_dt_*.mat`, `EP_LGIM_dt_*.mat`).

   Note: `.mat` outputs are not tracked in git (see repository `.gitignore`).

   Computational time: each integrator saves `executionTime` (integration only; see main
   `README.md`). At the top of `Run_Pendulum_Batch.m`, `n_timing_repeats` controls how many
   times each integration is run; the median is stored. Default is **1**. Use **5** for
   smoother timing curves; expect integration cost to scale by about **×5** (full pendulum
   batch with 5 repeats was on the order of ~1 h on a typical laptop in testing).

   Estimated batch runtime (default `n_timing_repeats = 1`): **~10 minutes** on the system
   described above (plus a few minutes for `PendulumPlotting.m`).

4) After the batch has finished, run `PendulumPlotting.m` to compare methods against
   `ExactRef_HighAccuracy.mat` and generate figures.

For a faster test run, reduce `time_steps_ms` inside `Run_Pendulum_Batch.m` (for example
`time_steps_ms = [10, 50, 100];` instead of the full list), similar to the example in
`READMEHeavyTop.txt`.
