function [gain, covariance, metadata] = static_fast_gain (deltaTRk, ~, ~, ~, ~, ~, ~)
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
    covariance = eye(11);

    metadata = [];
end
