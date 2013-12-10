function [gain, covariance] = kalman_time_varying_gain (~, A, G, C, Q, R, P)
    covariance_a_priori = A * P * A' + G * Q * G';
    gain = covariance_a_priori * C' / (C * covariance_a_priori * C' + R);
    covariance = (eye(11) - gain * C) * covariance_a_priori;
end