function state = create_state (x, time, innovation, covariance, Phi, residual, process_gain, gain, metadata)
    state = [];
    state.time = time;
    state.position = x(1:3);
    state.velocity = x(4:6);
    state.acceleration = x(7:9);
    state.clock_offset = x(10);
    state.clock_rate_offset = x(11);
    state.innovation = innovation;
    state.covariance = covariance;
    state.Phi = Phi;
    state.residual = residual;
    state.process_gain = process_gain;
    state.gain = gain;
    state.metadata = metadata;
    
    if length(x) > 11
        state.pseudorange_bias = x(12:43);
        state.doppler_shift_bias = x(44:75);
    end
end
