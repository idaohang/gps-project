function retval = sum_sym_kl_divergence(x_hat1, x_hat2)
    assert (length(x_hat1) == length(x_hat2));
    len = length(x_hat1);
    retval = zeros(len, 1);
    for i=1:len
        retval(i) = sym_kl_divergence(x_hat1{i}, x_hat2{i});
    end
    retval = sum(retval);
end