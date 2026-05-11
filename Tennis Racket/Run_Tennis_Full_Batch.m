% Run all Tennis Racket integrators with batch-controlled time steps (ms).
%
% Order: (1) ode45 reference -> TennisRef.mat
%        (2) each loop sets dt (seconds), save_filename, then run(...)
% Same convention as Heavy Top Run_HeavyTop_batch.m (reference + batch).

%% Reference solution (ode45) — TennisRef.mat
run('Tennis_ode45');

%% Full range (ms) — EP LGIM, NE EP, wrapped Cartesian NE, Cartesian LGIM
time_steps_ms = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0];
%time_steps_ms = [100.0, 200.0, 500.0];

% EP LGIM (Tennis_EP_LGIM.m)
for j = 1:numel(time_steps_ms)
    dt = time_steps_ms(j) / 1000;
    save_filename = sprintf('EP_LGIM_dt_%0.1fms.mat', time_steps_ms(j));
    run('Tennis_EP_LGIM');
end

% Newton–Euler EP (Tennis_NE_EP.m)
for j = 1:numel(time_steps_ms)
    dt = time_steps_ms(j) / 1000;
    save_filename = sprintf('EP_NE_dt_%0.1fms.mat', time_steps_ms(j));
    run('Tennis_NE_EP');
end

% Wrapped Cartesian NE (wrap_Tennis_NE_Cartesian.m)
for j = 1:numel(time_steps_ms)
    dt = time_steps_ms(j) / 1000;
    save_filename = sprintf('wrCart_NE_dt_%0.1fms.mat', time_steps_ms(j));
    run('wrap_Tennis_NE_Cartesian');
end

% Cartesian LGIM (Tennis_Cart_LGIM.m)
for j = 1:numel(time_steps_ms)
    dt = time_steps_ms(j) / 1000;
    save_filename = sprintf('wrLGIM_dt_%0.1fms.mat', time_steps_ms(j));
    run('Tennis_Cart_LGIM');
end

%% Classic EP (Tennis_Classic_EP.m): 0.1 ms … 20 ms only
time_steps_ms_classic_ep = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0];

for j = 1:numel(time_steps_ms_classic_ep)
    dt = time_steps_ms_classic_ep(j) / 1000;
    save_filename = sprintf('ClassicEP_dt_%0.1fms.mat', time_steps_ms_classic_ep(j));
    run('Tennis_Classic_EP');
end

%% Classic Cartesian (Tennis_Classic_Cart.m): 0.1 ms and 0.125 ms only
time_steps_ms_classic_cart = [0.1, 0.125];

for j = 1:numel(time_steps_ms_classic_cart)
    dt = time_steps_ms_classic_cart(j) / 1000;
    save_filename = sprintf('ClassicCart_dt_%gms.mat', time_steps_ms_classic_cart(j));
    run('Tennis_Classic_Cart');
end
