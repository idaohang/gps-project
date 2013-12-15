function ts = nis (x_hat)
    x_hat = x_hat(2:end);
    ts = zeros(size(x_hat));
    
    for i=1:size(x_hat)
        v = x_hat{i}.innovation;
        S = x_hat{i}.metadata.S;
        ts(i) = v(1:8)' / S(1:8,1:8) * v(1:8);
    end
    
    ts = ts(ts ~= 0);
end
