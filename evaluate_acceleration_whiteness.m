function evaluate_acceleration_whiteness(x_hat)
    accelerations = zeros(3, length(x_hat));
    for i=1:length(x_hat)
        accelerations(:,i) = x_hat{i}.acceleration;
    end
    [h, p_value, R] = whiteness_test (accelerations, 100);
    
end