function retval = kl_divergence (state1, state2)
    diff = state1.position - state2.position;
    cov1 = state1.covariance(1:3,1:3);
    cov2 = state2.covariance(1:3,1:3);
    retval = 0.5 * (log(det(cov1)/det(cov2)) + trace(cov1 \ cov2) + diff' / cov1 * diff - 3);
end