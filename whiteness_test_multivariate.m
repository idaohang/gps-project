function [h, p_value, Q] = whiteness_test_multivariate (x, maxlag, alpha)
    N = size(x, 2);
    
    if nargin < 2
        maxlag = min(120, floor(N/2)-2);
    end
    
    if nargin < 3
        alpha = 0.05;
    end
    
    K = maxlag;
    
    %x = bsxfun(@minus, x, mean(x, 2));
    
    Y0 = zeros(size(x,1));
    for i=1:N
        Y0 = Y0 + x(:,i) * x(:,i)';
    end
    Y0 = Y0/N;

    Y = cell(K, 1);
    for k=1:K
        Y{k} = zeros(size(x,1));
        
        for i=1:N-k
            Y{k} = Y{k} + x(:,i) * x(:,i+k)';
        end
        Y{k} = Y{k}/sqrt(N*(N-k));
    end

    Gamma = cell(K, 1);
    for k=1:K
        Gamma{k} = zeros(size(x,1));

        for i=1:size(Gamma{k},1)
            for j=1:size(Gamma{k},2)
                Gamma{k}(i,j) = Y{k}(i,j)/sqrt(Y0(i,i)*Y0(j,j));
            end
        end
    end

    Psi = zeros(size(x,1), 1);
    for j=1:size(Psi)
        for i=1:K
            Psi(j) = Psi(j) + Gamma{i}(j,j)^2;
        end
    end

    test_statistics = N * Psi;
    Q = chi2inv(1-alpha, K);
    p_value = 1-chi2cdf(test_statistics, K);
    h = test_statistics > Q;
    disp(test_statistics')
    disp(h')
end
