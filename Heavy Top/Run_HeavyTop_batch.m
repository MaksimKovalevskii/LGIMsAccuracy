% Single batch: run all Heavy Top programs for a list of time steps (in ms),
% plus a single reference solution with ode45 (fixed time step).
%
% Time steps in milliseconds:
%time_steps_ms = [3, 4, 10];
time_steps_ms = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 1.5, 2, 3, 4, 10];

% EP Lie group integrator
for j = 1:numel(time_steps_ms)
    dt = time_steps_ms(j) / 1000; % convert ms -> seconds
    save_filename = sprintf('HT_EPLGIM_dt_%.2fms.mat', time_steps_ms(j));
    run('HT_EP_LGIM');
end

% EP Newton-Euler
for j = 1:numel(time_steps_ms)
    dt = time_steps_ms(j) / 1000; % convert ms -> seconds
    save_filename = sprintf('HTNE_EP_dt_%.2fms.mat', time_steps_ms(j));
    run('HT_EP_NE');
end

% Cartesian Newton-Euler
for j = 1:numel(time_steps_ms)
    dt = time_steps_ms(j) / 1000; % convert ms -> seconds
    save_filename = sprintf('HT_CartNE_dt_%.2fms.mat', time_steps_ms(j));
    run('WrapHTCartNE');
end

% Cartesian Lie group integrator
for j = 1:numel(time_steps_ms)
    dt = time_steps_ms(j) / 1000; % convert ms -> seconds
    save_filename = sprintf('HT_CartLGIM_dt_%.2fms.mat', time_steps_ms(j));
    run('NewWrapHTCartLGIM');
end

% Classic EP (second-order ODE formulation)
for j = 1:numel(time_steps_ms)
    dt = time_steps_ms(j) / 1000; % convert ms -> seconds
    save_filename = sprintf('HTClassicEP_dt_%.2fms.mat', time_steps_ms(j));
    run('HT_ClassicEP');
end

% % Reference solution with ode45 (fixed, fine time step inside HT_ode45)
% save_filename = 'HTRef.mat';
% run('HT_ode45');
