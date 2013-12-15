function Q = process_covariance (process_gain, bias)
    if bias
        Q = process_gain * eye(68);
        Q(5:36,5:36) = 10*process_gain*eye(32);
        Q(37:68,37:68) = process_gain*eye(32);
    else
        Q= process_gain * eye(4);
    end
end