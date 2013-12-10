function retval = sym_kl_divergence (state1, state2)
    retval = kl_divergence(state1, state2) + kl_divergence(state2, state1);
end