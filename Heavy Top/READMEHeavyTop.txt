READMEHeavyTop

The heavy top is a classical benchmark problem in multibody system dynamics and a standard
test case for Lie group integration methods.

System used for development and testing:
- Operating system: Windows 10 (64‑bit)
- Machine type: Work laptop, multi‑core CPU
- MATLAB: R2024b (64‑bit), standard desktop installation
- Disk: At least 25 GB free space recommended for full batch

Instructions

1) Copy all files from this `Heavy Top` directory into a single folder on your local machine (or clone the GitHub repository and work directly in the `Heavy Top` folder).
2) Open MATLAB R2024b (or a compatible recent version) and set the current folder to this `Heavy Top` directory.
3) Run `Run_HeavyTop_batch.m`.  
   This script will:
   - Simulate the heavy top problem with 5 numerical integration methods, for a list of time steps given in milliseconds.
   - For each method and each time step, save a `.mat` database with simulation results.
   - Additionally compute a high‑accuracy reference solution using `ode45` and save it as `HTRef.mat`.
   Note: running the full batch requires at least 25 GB of free disk space.

Estimated runtime for full batch (all methods, all time steps): 30–45 minutes on a typical modern laptop similar to the system described above. Actual time depends on CPU performance and MATLAB version.

For a faster test run, use only larger time steps:
- In `Run_HeavyTop_batch.m`, set for example  
  `time_steps_ms = [3, 4, 10];`  
  instead of  
  `time_steps_ms = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 1.5, 2, 3, 4, 10];`

4) After the batch has finished, run `HeavyTop_plotting.m` to compare methods and generate plots.  
   Estimated runtime for the plotting script: ~10 minutes on the system described above.