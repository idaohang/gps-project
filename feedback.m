function next_state = feedback (deltaTRk, state, measurements, gain_mode)
    R = measurement_covariance (measurements);
    
    Q = process_covariance ();
    
    A = process_matrix (deltaTRk);
    
    G = zeros(11, 4);
    G(7:9,1:3) = eye(3);
    G(11, 4) = 1;
    
    C = zeros(8, 11);
    C(1:6, 1:6) = eye(6);
    C(7:8, 10:11) = eye(2);
    
    z = measurement_vector (measurements);
    x = state_vector (state);

    gain_fcn = gain_factory (gain_mode);
    
    [gain, covariance] = gain_fcn (deltaTRk, A, G, C, Q, R, state.covariance);
    
    prediction = A * x;
    innovation = z - C * prediction;
    x = prediction + gain * innovation;
    
    Phi = (eye(11) - gain * C) * A;
    
    next_state = create_state (x, measurements.time, innovation, covariance, Phi);
end