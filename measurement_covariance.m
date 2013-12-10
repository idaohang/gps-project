function R = measurement_covariance (measurements)
    constant;
    
    R = zeros(8, 8);
    R(1:3, 1:3) = measurements.sigmaPR^2*measurements.Q(1:3, 1:3);
    R(1:3, 7) = measurements.sigmaPR^2*measurements.Q(1:3, 4);
    R(7, 1:3) = measurements.sigmaPR^2*measurements.Q(4, 1:3);
    R(7, 7) = measurements.sigmaPR^2*measurements.Q(4, 4);
    R(4:6, 4:6) = (lambdaL1*measurements.sigmaDopp)^2*measurements.Qv(1:3, 1:3);
    R(4:6, 8) = (lambdaL1*measurements.sigmaDopp)^2*measurements.Qv(1:3, 4);
    R(8, 4:6) = (lambdaL1*measurements.sigmaDopp)^2*measurements.Qv(4, 1:3);
    R(8, 8) = (lambdaL1*measurements.sigmaDopp)^2*measurements.Qv(4, 4);
end