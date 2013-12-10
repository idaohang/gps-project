function gain_fcn = gain_factory (gain_mode)
    if strcmp(gain_mode, 'static-slow')
        gain_fcn = @static_slow_gain;
    elseif strcmp(gain_mode, 'static-fast')
        gain_fcn = @static_fast_gain;
    elseif strcmp(gain_mode, 'kalman-steady-state')
        gain_fcn = @kalman_steady_state_gain;
    elseif strcmp(gain_mode, 'kalman-time-varying')
        gain_fcn = @kalman_time_varying_gain;
    end
end