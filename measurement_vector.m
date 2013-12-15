function z = measurement_vector (measurements)
    z = [
        measurements.position; 
        measurements.velocity; 
        measurements.clock_offset; 
        measurements.clock_rate_offset
    ];

    if isfield(measurements, 'pseudorange_bias')
        z = [z; measurements.pseudorange_bias(:,2); measurements.doppler_shift_bias(:,2)];
    end
end