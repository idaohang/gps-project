whiteness_test.mfunction evaluate_jerk_whiteness(x_hat)
    jerk = compute_jerk (x_hat);
    [h, p_value, R] = whiteness_test (jerk);
    
end