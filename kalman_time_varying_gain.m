function [gain, covariance, metadata] = kalman_time_varying_gain (~, A, G, C, Q, R, P)
    covariance_a_priori = A * P * A' + G * Q * G';
    S = C * covariance_a_priori * C' + R;
    gain = covariance_a_priori * C' / S;
    covariance = (eye(size(gain,1)) - gain * C) * covariance_a_priori;
    metadata.S = S;
end
