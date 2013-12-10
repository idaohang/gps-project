function next_state = static_feedback(deltaTRk, state, measurements, ...
                                      gain_mode)
    prediction = [
        state.position + ...
        deltaTRk * state.velocity + ...
        0.5 * deltaTRk^2 * state.acceleration;
        
        state.velocity + ...
        deltaTRk * state.acceleration;
        
        state.acceleration;
        
        state.clock_offset + ...
        deltaTRk * state.clock_rate_offset;
        
        state.clock_rate_offset;
    ];

    z = [
        measurements.position;
        measurements.velocity;
        measurements.clock_offset;
        measurements.clock_rate_offset;
    ];

    innovation = z - prediction([1:6 10:11]);

    if strcmp(gain_mode, 'slow')
        gain = zeros(11, 8);
        gain(1:3,1:3) = 8.44567e-3 * eye(3);
        gain(1:3,4:6) = 0.508214 * deltaTRk * eye(3);
        
        gain(4:6,1:3) = 3.67184e-5/deltaTRk * eye(3);
        gain(4:6,4:6) = 0.975104 * eye(3);
        
        gain(7:9,1:3) = -3.29082e-5/deltaTRk^2 * eye(3);
        gain(7:9,4:6) = 0.927802/deltaTRk * eye(3);
        
        gain(10, 7) = 8.462e-3;
        gain(10, 8) = 0.0271035 * deltaTRk;
        
        gain(11, 7) = 1.95823e-6/deltaTRk;
        gain(11, 8) = 0.972659;
    elseif strcmp(gain_mode, 'fast')
        gain = zeros(11, 8);
        gain(1:3,1:3) = 8.44612e-3 * eye(3);
        gain(1:3,4:6) = 0.495778 * deltaTRk * eye(3);
        
        gain(4:6,1:3) = 3.582e-5/deltaTRk * eye(3);
        gain(4:6,4:6) = 0.999979 * eye(3);
        
        gain(7:9,1:3) = -3.55164e-5/deltaTRk^2 * eye(3);
        gain(7:9,4:6) = 1.000009/deltaTRk * eye(3);
        
        gain(10, 7) = 8.46395e-3;
        gain(10, 8) = 2.86552e-6 * deltaTRk;
        
        gain(11, 7) = 2.07034e-10/deltaTRk;
        gain(11, 8) = 0.999997;
    end
    
    x = prediction + gain * innovation;

    next_state = [];
    next_state.innovation = innovation;
    next_state.time = measurements.time;
    next_state.position = x(1:3);
    next_state.velocity = x(4:6);
    next_state.acceleration = x(7:9);
    next_state.clock_offset = x(10);
    next_state.clock_rate_offset = x(11);
end