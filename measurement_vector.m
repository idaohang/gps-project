function z = measurement_vector (measurements)
    z = [
        measurements.position; 
        measurements.velocity; 
        measurements.clock_offset; 
        measurements.clock_rate_offset
    ];
end