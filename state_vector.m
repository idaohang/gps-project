function x = state_vector (state)
    x = [
        state.position; 
        state.velocity; 
        state.acceleration; 
        state.clock_offset; 
        state.clock_rate_offset
    ];

    if isfield(state, 'pseudorange_bias')
        x = [x; state.pseudorange_bias; state.doppler_shift_bias];
    end
end