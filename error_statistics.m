function [m, s] = error_statistics (estimated, correct)
    assert(length(estimated) == length(correct));
    len = length(estimated);
    error = zeros(3, len);
    for i=1:len
        error(:, i) = estimated{i}.position - correct{i}.position;
    end
    convert to VEN
    m = mean(error);
    s = sqrt(var(error));
end