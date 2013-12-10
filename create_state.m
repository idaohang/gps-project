function state = create_state (x, time, innovation, covariance, Phi)
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
end