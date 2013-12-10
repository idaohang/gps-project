function plot_max_abs_eigenvalues (x_hat)
    eigenvalues = zeros(length(x_hat)-1, 1);
    for i=2:length(x_hat)
        e = eig(x_hat{i}.Phi);
        eigenvalues(i-1) = max(abs(e));
    end
    scatter(1:length(eigenvalues), eigenvalues, 'rx')
end