function x = state_vector (state)
    x = [
        state.position; 
        state.velocity; 
        state.acceleration; 
        state.clock_offset; 
        state.clock_rate_offset
    ];
end