function R = measurement_covariance (measurements, bias)
    constant;

    if bias
        num_measurements = 8+size(measurements.pseudorange_bias,1)+size(measurements.doppler_shift_bias,1);
    else
        num_measurements = 8;
    end

    R = zeros(num_measurements);
    R(1:3, 1:3) = measurements.sigmaPR^2*measurements.Q(1:3, 1:3);
    R(1:3, 7) = measurements.sigmaPR^2*measurements.Q(1:3, 4);
    R(7, 1:3) = measurements.sigmaPR^2*measurements.Q(4, 1:3);
    R(7, 7) = measurements.sigmaPR^2*measurements.Q(4, 4);
    R(4:6, 4:6) = (lambdaL1*measurements.sigmaDopp)^2*measurements.Qv(1:3, 1:3);
    R(4:6, 8) = (lambdaL1*measurements.sigmaDopp)^2*measurements.Qv(1:3, 4);
    R(8, 4:6) = (lambdaL1*measurements.sigmaDopp)^2*measurements.Qv(4, 1:3);
    R(8, 8) = (lambdaL1*measurements.sigmaDopp)^2*measurements.Qv(4, 4);
    
    if bias
        R(9:8+size(measurements.pseudorange_bias,1),9:8+size(measurements.pseudorange_bias,1)) = diag(measurements.pseudorange_var(:,2));
        R(9+size(measurements.pseudorange_bias,1):end,9+size(measurements.pseudorange_bias,1):end) = diag(measurements.doppler_shift_var(:,2));
    end
end
