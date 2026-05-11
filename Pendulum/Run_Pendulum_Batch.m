function Run_Pendulum_Batch()
%RUN_PENDULUM_BATCH  Reference solution + six integrators at shared time steps.
%
% Writes .mat files compatible with PendulumPlotting.m naming:
%   ExactRef_HighAccuracy.mat, EP_dt_*.mat, Cart_dt_*.mat, EP_NE_dt_*.mat,
%   wrCart_NE_dt_*.mat, wrLGIM_dt_*.mat, EP_LGIM_dt_*.mat

    here = fileparts(mfilename('fullpath'));
    oldPwd = cd(here);
    try
        time_steps_ms = [0.2, 0.5, 1, 2, 5, 10, 20, 50, 100];

        % --- Reference (high-accuracy ode45) ---
        clear save_filename
        run('ThetaRefode45');

        % Order aligned with PendulumPlotting file_patterns / typical workflow
        integrators = {
            'ClassicEP100',             'EP_dt_%.1fms.mat'
            'wrap_ClassicCartesian100', 'Cart_dt_%.1fms.mat'
            'EP_NE_100',                'EP_NE_dt_%.1fms.mat'
            'wrapCartNE100',            'wrCart_NE_dt_%.1fms.mat'
            'wrap_Cart_LGIM_100',       'wrLGIM_dt_%.1fms.mat'
            'EP_LGIM_100',              'EP_LGIM_dt_%.1fms.mat'
            };

        for j = 1:numel(time_steps_ms)
            dt_ms = time_steps_ms(j);
            dt = dt_ms / 1000;

            for k = 1:size(integrators, 1)
                scriptName = integrators{k, 1};
                pat = integrators{k, 2};
                save_filename = sprintf(pat, dt_ms);
                run(scriptName);
            end
        end
    finally
        cd(oldPwd);
    end
end
