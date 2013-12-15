function [gain, covariance,t] = kalman_steady_state_gain (~, A, G, C, Q, R, ~,~)
    [gain, ~, covariance] = dlqe(A, G, C, Q, R);
    t = 0;
end