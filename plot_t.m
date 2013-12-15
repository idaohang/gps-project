function pass = plot_t (x_hat)
    x_hat = x_hat(2:end);
    ts = zeros(size(x_hat));
    
    for i=1:size(x_hat)
        v = x_hat{i}.innovation;
        S = x_hat{i}.metadata.S;
        ts(i) = v(1:8)' / S(1:8,1:8) * v(1:8);
    end
    
    ts = ts(ts ~= 0);
    
    %ts = ts(1:4:end);
    semilogy(1:length(ts), ts, 'rx')
    
    t = sum(ts)/length(ts);
    
    n_v = 8;
    N = length(ts);
    
    alpha = 0.01;
    r1 = chi2inv( alpha/2,  N*n_v ) / N;
    r2 = chi2inv( 1 - alpha/2,  N*n_v ) / N;
    pass = (t > r1) && (t < r2);
    disp([r1 t r2])
end
