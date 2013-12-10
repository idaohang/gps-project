function [gain, covariance] = static_slow_gain (deltaTRk, ~, ~, ~, ~, ~, ~)
    gain = zeros(11, 8);
    gain(1:3,1:3) = 8.44567e-3 * eye(3);
    gain(1:3,4:6) = 0.508214 * deltaTRk * eye(3);

    gain(4:6,1:3) = 3.67184e-5/deltaTRk * eye(3);
    gain(4:6,4:6) = 0.975104 * eye(3);

    gain(7:9,1:3) = -3.29082e-5/deltaTRk^2 * eye(3);
    gain(7:9,4:6) = 0.927802/deltaTRk * eye(3);

    gain(10, 7) = 8.462e-3;
    gain(10, 8) = 0.0271035 * deltaTRk;

    gain(11, 7) = 1.95823e-6/deltaTRk;
    gain(11, 8) = 0.972659;
    covariance = eye(11);
end