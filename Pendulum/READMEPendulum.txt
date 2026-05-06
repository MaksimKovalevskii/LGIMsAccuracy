READMEPendulum

The pendulum example is a constrained multibody benchmark with smooth, large-amplitude rotation.
It is used here to compare six integration formulations (conventional and Newton–Euler style,
with Euler parameters and Cartesian rotation vectors, plus LGIM variants) against a
high-accuracy numerical reference.

System used for development and testing:
- Operating system: Windows 10/11 (64-bit)
- Machine type: Work laptop, multi-core CPU
- MATLAB: R2024b (64-bit), standard desktop installation
- Disk: At least 5 GB free space recommended for a full batch of all methods and time steps

Instructions

1) Copy all files from this `Pendulum` directory into a single folder on your local machine
   (or clone the GitHub repository and work directly in the `Pendulum` folder).
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

4) After the batch has finished, run `PendulumPlotting.m` to compare methods against
   `ExactRef_HighAccuracy.mat` and generate figures. Estimated runtime for the plotting script:
   on the order of minutes on the system described above (depends on disk I/O and figure count).

For a faster test run, reduce `time_steps_ms` inside `Run_Pendulum_Batch.m` (for example
`time_steps_ms = [10, 50, 100];` instead of the full list), similar to the example in
`READMEHeavyTop.txt`.
