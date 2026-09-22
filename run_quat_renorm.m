%% Re-run only Conv-Quat and NE-Quat (comment 12 renormalization)
% Leaves Cartesian / LGIM / ode45 .mat files untouched. Overwrites the
% existing Conv-Quat and NE-Quat databases, then rebuilds the comparison plots.
% Time-step lists match Run_Pendulum_Batch, Run_Tennis_Full_Batch, Run_HeavyTop_batch.

repo = fileparts(mfilename('fullpath'));
addpath(repo);
n_timing_repeats = 1;
do_plotting = true;

%% Pendulum
disp_time("Pendulum: Conv-Quat + NE-Quat")
clearvars -except repo n_timing_repeats do_plotting
clear global
cd(fullfile(repo, "Pendulum"))
time_steps_ms = [0.2, 0.5, 1, 2, 5, 10, 20, 50, 100];
integrators = {
    'ClassicEP100', 'EP_dt_%.1fms.mat'
    'EP_NE_100',    'EP_NE_dt_%.1fms.mat'
    };
for j = 1:numel(time_steps_ms)
    dt_ms = time_steps_ms(j);
    dt = dt_ms / 1000;
    for k = 1:size(integrators, 1)
        save_filename = sprintf(integrators{k, 2}, dt_ms);
        run(integrators{k, 1});
    end
end
if do_plotting
    PendulumPlotting;
end
cd(repo)

%% Tennis racket
disp_time("Tennis racket: Conv-Quat + NE-Quat")
clearvars -except repo n_timing_repeats do_plotting
clear global
cd(fullfile(repo, "Tennis Racket"))

time_steps_ms = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0];
for j = 1:numel(time_steps_ms)
    dt = time_steps_ms(j) / 1000;
    save_filename = sprintf('EP_NE_dt_%0.1fms.mat', time_steps_ms(j));
    run('Tennis_NE_EP');
end

time_steps_ms_classic_ep = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0];
for j = 1:numel(time_steps_ms_classic_ep)
    dt = time_steps_ms_classic_ep(j) / 1000;
    save_filename = sprintf('ClassicEP_dt_%0.1fms.mat', time_steps_ms_classic_ep(j));
    run('Tennis_Classic_EP');
end
if do_plotting
    Tennis_Plotting;
end
cd(repo)

%% Heavy top
disp_time("Heavy top: Conv-Quat + NE-Quat")
clearvars -except repo n_timing_repeats do_plotting
clear global
cd(fullfile(repo, "Heavy Top"))
time_steps_ms = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 1.5, 2, 3, 4, 10];
for j = 1:numel(time_steps_ms)
    dt = time_steps_ms(j) / 1000;
    save_filename = sprintf('HTNE_EP_dt_%.2fms.mat', time_steps_ms(j));
    run('HT_EP_NE');
end
for j = 1:numel(time_steps_ms)
    dt = time_steps_ms(j) / 1000;
    save_filename = sprintf('HTClassicEP_dt_%.2fms.mat', time_steps_ms(j));
    run('HT_ClassicEP');
end
if do_plotting
    HeavyTop_plotting;
end
cd(repo)

disp_time("Quat renormalization runs end")

function disp_time(X)
    disp(X);
    disp(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
end
