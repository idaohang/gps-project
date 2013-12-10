function jerk = compute_jerk (x_hat)
    jerk = zeros(3, length(x_hat)-1);
    for i=1:length(x_hat)-1
        jerk(:,i) = (x_hat{i+1}.acceleration - x_hat{i}.acceleration)/(x_hat{i+1}.time - x_hat{i}.time);
    end
end