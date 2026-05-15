%% Run Pendulum Example - reference + plots
disp_time("Run Pendulum Example - reference + plots")
clear;
cd("Pendulum/")
Run_Pendulum_Batch;
PendulumPlotting;
cd("..")

%% Tennis racket
disp_time("Run Tennis Racket - reference + plots")
clear;
cd("Tennis Racket/")
Run_Tennis_Full_Batch;
Tennis_Plotting;
cd("..")

%% Heavy Top Example
disp_time("Run Heavy Top - reference + plots")
clear;
cd("Heavy Top/")
Run_HeavyTop_batch;
HeavyTop_plotting;
cd("..")

%% End the whole run
disp_time("All runs end")


%% Helper to see time and not infer with internal tic/toc
function disp_time(X)
    disp(X);
    disp(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
end