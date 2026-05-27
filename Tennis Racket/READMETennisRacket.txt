READMETennisRacket

The tennis racket theorem, also known as the intermediate axis theorem or Dzhanibekov effect, describes the unstable rotational behavior of a rigid body about its intermediate principal axis. When a tennis racket (or any object with three distinct principal moments of inertia) is thrown
with initial rotation predominantly about the intermediate axis, the body undergoes characteristic tumbling motion with periodic half-rotations, during which the object’s orientation flips through π radians.

System used for development and testing:
- Operating system: Windows 10/11 (64-bit)
- Machine type: Work laptop, multi-core CPU
- MATLAB: R2024b (64-bit), standard desktop installation
- Disk: Full batch `.mat` outputs are typically ~3 GB; allow **at least 4 GB** free (headroom).

Instructions

1) Copy all files from this `Tennis Racket` directory into a single folder on your local machine
   (or clone the GitHub repository and work directly in the `Tennis Racket` folder).
   Alternatively, from the repository root, run `run_all` to execute all three examples in one go
   (see the main `README.md`).
2) Open MATLAB R2024b (or a compatible recent version) and set the current folder to this
   `Tennis Racket` directory.
3) Run `Run_Tennis_Full_Batch.m`.
   This script will:
   - First run `Tennis_ode45.m` to build a high-accuracy reference and save it as `TennisRef.mat`.
   - Then simulate six formulations (Classic EP, Classic Cartesian, NE EP, wrapped NE Cartesian,
     Cartesian LGIM, EP LGIM) over the time-step lists defined in the batch script:
       * Full shared list (ms) for EP LGIM, NE EP, wrapped Cartesian NE, and Cartesian LGIM.
       * A shorter list (0.1–20 ms) for Classic EP.
       * Two steps (0.1 ms and 0.125 ms) for Classic Cartesian.
   - For each run, write a `.mat` file with the naming pattern used in `Tennis_Plotting.m`.

   Note: `.mat` outputs are not tracked in git (see repository `.gitignore`).

   Computational time: each integrator saves `executionTime` (integration only; see main
   `README.md`). At the top of `Run_Tennis_Full_Batch.m`, `n_timing_repeats` controls how many
   times each integration is run; the median is stored. Default is **1**. Use **5** for
   smoother timing curves; expect integration cost to scale by about **×5**.

   Estimated batch runtime (default `n_timing_repeats = 1`): **~10 minutes** on the system
   described above (plus a few minutes for `Tennis_Plotting.m`).

4) After the batch has finished, run `Tennis_Plotting.m` to compare methods against `TennisRef.mat`
   and generate figures.

For a faster test run, reduce the time-step vectors inside `Run_Tennis_Full_Batch.m` (comment the
full lists and use a small subset, similar to the Heavy Top example in `READMEHeavyTop.txt`).
