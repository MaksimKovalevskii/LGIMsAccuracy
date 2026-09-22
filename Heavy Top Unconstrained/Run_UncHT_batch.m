% Unconstrained heavy top: 5 methods + ode45 reference.
% Classic Cart excluded (singular at practical steps, same as tennis).
% Run from this folder. Overwrites .mat in this folder only.
n_timing_repeats = 1;

%time_steps_ms = [3, 4, 10];
time_steps_ms = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 1.5, 2, 3, 4, 10];

% EP Lie group integrator
for j = 1:numel(time_steps_ms)
    dt = time_steps_ms(j) / 1000;
    save_filename = sprintf('UncHT_EPLGIM_dt_%.2fms.mat', time_steps_ms(j));
    run('UncHT_EP_LGIM');
end

% EP Newton-Euler
for j = 1:numel(time_steps_ms)
    dt = time_steps_ms(j) / 1000;
    save_filename = sprintf('UncHT_NE_EP_dt_%.2fms.mat', time_steps_ms(j));
    run('UncHT_EP_NE');
end

% Cartesian Newton-Euler
for j = 1:numel(time_steps_ms)
    dt = time_steps_ms(j) / 1000;
    save_filename = sprintf('UncHT_CartNE_dt_%.2fms.mat', time_steps_ms(j));
    run('UncHT_CartNE');
end

% Cartesian Lie group integrator
for j = 1:numel(time_steps_ms)
    dt = time_steps_ms(j) / 1000;
    save_filename = sprintf('UncHT_CartLGIM_dt_%.2fms.mat', time_steps_ms(j));
    run('UncHT_CartLGIM');
end

% Classic EP (second-order ODE)
for j = 1:numel(time_steps_ms)
    dt = time_steps_ms(j) / 1000;
    save_filename = sprintf('UncHT_ClassicEP_dt_%.2fms.mat', time_steps_ms(j));
    run('UncHT_ClassicEP');
end

% Reference (unconstrained ode45)
save_filename = 'UncHTRef.mat';
run('UncHT_ode45');
