function R = process_covariance (process_gain)
    R = process_gain * eye(4);
end