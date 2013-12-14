function retval = evaluate_residuals_whiteness(x_hat)
    x_hat = x_hat(2:end);
    residuals = zeros(11, length(x_hat));
    lla = latlong(x_hat{1}.position');
    phi = lla(1)*pi/180;
    lambda = lla(2)*pi/180;
    A_VEN_ECEF = ...
      [cos(phi),0,sin(phi);0,1,0;-sin(phi),0,cos(phi)]*...
      [cos(lambda),sin(lambda),0;...
        -sin(lambda),cos(lambda),0;0,0,1];

    for i=1:length(x_hat)
        residuals(:,i) = x_hat{i}.residual(1:11);
    end
    residuals(1:3,:) = A_VEN_ECEF * residuals(1:3,:);
    residuals(4:6,:) = A_VEN_ECEF * residuals(4:6,:);
    h = whiteness_test_portmanteau (residuals(1:3,:));
    return

    pass = zeros(4,1);
    fail = zeros(4,1);
    n = size(residuals,2)-120;
    skip = max(1,floor(n/100));
    for i=1:skip:n
        [h, p_value] = whiteness_test_multivariate (residuals(:,i:i+120));
        if max(h([1:3])) == 0
            pass(1) = pass(1) + 1;
        else
            fail(1) = fail(1) + 1;
        end
        if max(h([4:6])) == 0
            pass(2) = pass(2) + 1;
        else
            fail(2) = fail(2) + 1;
        end
        if h(7) == 0
            pass(3) = pass(3) + 1;
        else
            fail(3) = fail(3) + 1;
        end
        if h(8) == 0
            pass(4) = pass(4) + 1;
        else
            fail(4) = fail(4) + 1;
        end
        disp((pass./(pass+fail))')
    end
    retval = pass./(pass+fail);
end
