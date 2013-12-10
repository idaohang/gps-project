function evaluate_innovation_whiteness(x_hat)
    innovations = zeros(8, length(x_hat));
    for i=1:length(x_hat)
        innovations(:,i) = x_hat{i}.innovation;
    end
    [h, p_value, R] = whiteness_test (innovations, 5000);
    
end