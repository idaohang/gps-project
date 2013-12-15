function reject = whiteness_test_portnameteau (x, h, alpha)
    T = size(x, 2);
    K = size(x, 1);
    p = 0;
    
    h = floor(T/2);

    C0 = zeros(K);
    for t=1:T
        C0 = C0 + x(:,t)*x(:,t)';
    end
    C0 = C0/T;

    C = cell(h, 1);
    for i=1:h
        C{i} = zeros(K);
        for t=i+1:T
            C{i} = C{i} + x(:,t)*x(:,t-i)';
        end
        C{i} = C{i}/T;
    end

    Qh = 0;
    for j=1:h
        Qh = Qh + trace(C{j}'*inv(C0)*C{j}*inv(C0))/(T-j);
    end
    Qh = Qh*T*(T+2);

    dof = K*K*(h-p);

    pp=1-chi2cdf(Qh,dof);

    reject = pp < 0.05;

    disp(Qh)
    disp(pp)
end
