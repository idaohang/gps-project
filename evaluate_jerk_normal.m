function h = evaluate_jerk_normal (xhat)
    jerk = compute_jerk(xhat);
    
    C = cov(jerk');
    m = mean(jerk, 2);
    Cb = bsxfun(@minus, jerk, m);
    
    n = size(jerk, 2);
    
    A = 0;
    for i=1:n
        for j=1:n
            A = A + (Cb(:,i)'/C*Cb(:,j))^3;
        end
    end
    A = A / (6*n);
    
    k = size(jerk, 1);
    
    B = 0;
    for i=1:n
        B = B + (Cb(:,i)'/C*Cb(:,i))^2;
    end
    B = B / n - k * (k + 2);
    B = sqrt(n / (8*k*(k+2))) * B;
    
   
    df = 1/6*k*(k+1)*(k+2);
    
    
    p1 = 1 - chi2cdf(A, df);
    
    p2 = 1 - normcdf(abs(B));
    
    p1 > 0.05 && p2 > 0.05;
end