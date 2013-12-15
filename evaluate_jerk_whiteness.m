function h = evaluate_jerk_whiteness(x_hat)
    jerk = compute_jerk (x_hat);
    h = whiteness_test_portmanteau(jerk, 10, 0.05);
end
