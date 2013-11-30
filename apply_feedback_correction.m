function x_kp1_hat = apply_feedback_correction (x_kp1_bar, v_kp1, delta_T_Rk_clock_rate_offset, Q, Qv, y_n, x_nm1, gain_mode)
    if strcmp(gain_mode, 'slow')
        x_kp1_hat = apply_slow_feedback_correction (x_kp1_bar, v_kp1, delta_T_Rk_clock_rate_offset);
    elseif strcmp(gain_mode, 'fast')
        x_kp1_hat = apply_fast_feedback_correction (x_kp1_bar, v_kp1, delta_T_Rk_clock_rate_offset);
    elseif strcmp(gain_mode, 'kalman')
        x_kp1_hat = apply_kalman_feedback_correction (x_kp1_bar, v_kp1, delta_T_Rk_clock_rate_offset, Q, Qv, y_n, x_nm1);
    end
end