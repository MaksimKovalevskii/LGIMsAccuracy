READMETennisRacket

The tennis racket theorem, also known as the intermediate axis theorem or Dzhanibekov effect, describes the unstable rotational behavior of a rigid body about its intermediate principal axis. When a tennis racket (or any object with three distinct principal moments of inertia) is thrown
with initial rotation predominantly about the intermediate axis, the body undergoes characteristic tumbling motion with periodic half-rotations, during which the object’s orientation flips through π radians.

System used for development and testing:
- Operating system: Windows 10/11 (64-bit)
- Machine type: Work laptop, multi-core CPU
- MATLAB: R2024b (64-bit), standard desktop installation
- Disk: At least 10 GB free space recommended for a full batch of all methods and time steps

Instructions

1) Copy all files from this `Tennis Racket` directory into a single folder on your local machine
   (or clone the GitHub repository and work directly in the `Tennis Racket` folder).
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

4) After the batch has finished, run `Tennis_Plotting.m` to compare methods against `TennisRef.mat`
   and generate figures. Estimated runtime for the plotting script: on the order of minutes on
   the system described above (depends on disk I/O and figure count).

For a faster test run, reduce the time-step vectors inside `Run_Tennis_Full_Batch.m` (comment the
full lists and use a small subset, similar to the Heavy Top example in `READMEHeavyTop.txt`).
