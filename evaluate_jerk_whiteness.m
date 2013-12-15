function retval = evaluate_jerk_whiteness(x_hat)
    jerk = compute_jerk (x_hat);
    pass = 0;
    fail = 0;
    n = size(jerk,2)-120;
    skip = max(1,floor(n/100));
    for i=1:skip:n
        [h, p_value] = whiteness_test_multivariate (jerk(:,i:i+120));
        if max(h) == 0
            pass = pass + 1;
        else
            fail = fail + 1;
        end
    end
    retval = pass/(pass+fail);
end
