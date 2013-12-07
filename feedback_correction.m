function next_state = feedback_correction (deltaTRk, state, ...
                                           measurements, method, gain_mode)
    if strcmp(method, 'static')
        next_state = static_feedback (deltaTRk, state, measurements, gain_mode);
    elseif strcmp(method, 'kalman')
        next_state = kalman_feedback (deltaTRk, state, measurements, gain_mode);
    end
end