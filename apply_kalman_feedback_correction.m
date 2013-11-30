function x_kp1_hat = apply_kalman_feedback_correction (x_kp1_bar, v_kp1, delta_T_Rk_clock_rate_offset, Q, Qv, y, x_nm1)
    constant;
    
    sigma_pr = 30;
    sigma_D = 10;

    R = zeros(8, 8);
    R(1:3, 1:3) = sigma_pr^2*Q(1:3, 1:3);
    R(1:3, 7) = sigma_pr^2*Q(1:3, 4);
    R(7, 1:3) = sigma_pr^2*Q(4, 1:3);
    R(7, 7) = sigma_pr^2*Q(4, 4);
    R(4:6, 4:6) = (lambdaL1*sigma_D)^2*Qv(1:3, 1:3);
    R(4:6, 8) = (lambdaL1*sigma_D)^2*Qv(1:3, 4);
    R(8, 4:6) = (lambdaL1*sigma_D)^2*Qv(4, 1:3);
    R(8, 8) = (lambdaL1*sigma_D)^2*Qv(4, 4);
    
    Q = 0.1*eye(4);
    
    G = zeros(11, 4);
    G(7:9,1:3) = eye(3);
    G(11, 4) = 1;
    
    A = eye(11, 11);
    A(1:3, 4:6) = delta_T_Rk_clock_rate_offset * eye(3);
    A(1:3, 7:9) = delta_T_Rk_clock_rate_offset^2 * eye(3);
    A(4:6, 7:9) = delta_T_Rk_clock_rate_offset * eye(3);
    A(10, 11) = delta_T_Rk_clock_rate_offset;
    
    C = zeros(8, 11);
    C(1:6, 1:6) = eye(6);
    C(7:8, 10:11) = eye(2);
    
    M = dlqe(A, G, C, Q, R);
    
    y_n = [y.position; y.velocity; y.clock_offset; y.clock_rate_offset];
    x_n_nm1 = [x_nm1.position; x_nm1.velocity; x_nm1.acceleration; x_nm1.clock_offset; x_nm1.clock_rate_offset];
    
    x_n_n = x_n_nm1 + M * (y_n - C * x_n_nm1);
    x_np1_n = A * x_n_n;
    
    x_kp1_hat = [];
    x_kp1_hat.time = x_kp1_bar.time;
    x_kp1_hat.position = x_np1_n(1:3);
    x_kp1_hat.velocity = x_np1_n(4:6);
    x_kp1_hat.acceleration = x_np1_n(7:9);
    x_kp1_hat.clock_offset = x_np1_n(10);
    x_kp1_hat.clock_rate_offset = x_np1_n(11);
end