function retval = kl_divergence (state1, state2)
    x1 = state_vector (state1);
    x2 = state_vector (state2);
    diff = x1 - x2;
    cov1 = state1.covariance;
    cov2 = state2.covariance;
    retval = 0.5 * (log(det(cov1)/det(cov2)) + trace(cov1 \ cov2) + diff' / cov1 * diff - length(diff));
end