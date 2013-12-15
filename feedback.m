function next_state = feedback (deltaTRk, state, measurements, gain_mode, process_gain, measurement_gain)
    bias = isfield(state, 'pseudorange_bias');
    
    R = measurement_gain*measurement_covariance (measurements, bias);

    A = process_matrix (deltaTRk, bias);
    
    if bias
        num_states = 75;
        num_unobservable = 68;
    else
        num_states = 11;
        num_unobservable = 4;
    end
    num_measurements = size(R,1);
    G = zeros(num_states, num_unobservable);
    G(7:9,1:3) = eye(3);
    G(11, 4) = 1;
    if bias
        G(12:end,5:end) = eye(64);
    end
    
    C = zeros(num_measurements, num_states);
    C(1:6, 1:6) = eye(6);
    C(7:8, 10:11) = eye(2);
    if bias
        for i=1:size(measurements.pseudorange_bias,1)
            C(8+i, 11+measurements.pseudorange_bias(i,1)) = 1;
        end
        for i=1:size(measurements.doppler_shift_bias,1)
            C(8+size(measurements.pseudorange_bias,1)+i, 11+32+measurements.doppler_shift_bias(i,1)) = 1;
        end
    end

    z = measurement_vector (measurements);
    x = state_vector (state);

    gain_fcn = gain_factory (gain_mode);

    prediction = A * x;
    innovation = z - C * prediction;
    
    Q = process_covariance (process_gain, bias);
    
    [gain, covariance, metadata] = gain_fcn (deltaTRk, A, G, C, Q, R, state.covariance);
    
    x = prediction + gain * innovation;
    residual = prediction - x;
    
    Phi = (eye(num_states) - gain * C) * A;
    
    next_state = create_state (x, measurements.time, innovation, covariance, Phi, residual, process_gain, gain, metadata);

end
