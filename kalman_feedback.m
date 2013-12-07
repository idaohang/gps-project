function next_state = kalman_feedback (deltaTRk, state, ...
                                       measurements, gain_mode)
    constant;
    
    sigma_pr = 30;
    sigma_D = 10;

    R = zeros(8, 8);
    R(1:3, 1:3) = sigma_pr^2*measurements.Q(1:3, 1:3);
    R(1:3, 7) = sigma_pr^2*measurements.Q(1:3, 4);
    R(7, 1:3) = sigma_pr^2*measurements.Q(4, 1:3);
    R(7, 7) = sigma_pr^2*measurements.Q(4, 4);
    R(4:6, 4:6) = (lambdaL1*sigma_D)^2*measurements.Qv(1:3, 1:3);
    R(4:6, 8) = (lambdaL1*sigma_D)^2*measurements.Qv(1:3, 4);
    R(8, 4:6) = (lambdaL1*sigma_D)^2*measurements.Qv(4, 1:3);
    R(8, 8) = (lambdaL1*sigma_D)^2*measurements.Qv(4, 4);
    
    Q = 10*eye(4);
    
    G = zeros(11, 4);
    G(7:9,1:3) = eye(3);
    G(11, 4) = 1;
    
    A = eye(11, 11);
    A(1:3, 4:6) = deltaTRk * eye(3);
    A(1:3, 7:9) = deltaTRk^2 * eye(3);
    A(4:6, 7:9) = deltaTRk * eye(3);
    A(10, 11) = deltaTRk;
    
    C = zeros(8, 11);
    C(1:6, 1:6) = eye(6);
    C(7:8, 10:11) = eye(2);
    
    
    z = [measurements.position; measurements.velocity; measurements.clock_offset; measurements.clock_rate_offset];
    x = [state.position; state.velocity; state.acceleration; state.clock_offset; state.clock_rate_offset];

    if strcmp(gain_mode, 'steady-state')
        [kalman_gain, ~, covariance] = dlqe(A, G, C, Q, R);
    elseif strcmp(gain_mode, 'time-varying')
        covariance_a_priori = A * state.covariance * A' + G * Q * G';
        kalman_gain = covariance_a_priori * C' / (C * covariance_a_priori * C' + R);
        covariance = (eye(11) - kalman_gain * C) * covariance_a_priori;
    end
    
    prediction = A * x;
    innovation = z - C * prediction;
    x = prediction + kalman_gain * innovation;
    
    next_state = [];
    next_state.time = measurements.time;
    next_state.position = x(1:3);
    next_state.velocity = x(4:6);
    next_state.acceleration = x(7:9);
    next_state.clock_offset = x(10);
    next_state.clock_rate_offset = x(11);
    next_state.innovation = innovation;
    next_state.covariance = covariance;
end