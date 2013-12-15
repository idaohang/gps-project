function A = process_matrix (deltaTRk, bias)
    if bias
        A = eye(75);
    else
        A = eye(11);
    end
    A(1:3, 4:6) = deltaTRk * eye(3);
    A(1:3, 7:9) = 0.5 * deltaTRk^2 * eye(3);
    A(4:6, 7:9) = deltaTRk * eye(3);
    A(10, 11) = deltaTRk;
end
