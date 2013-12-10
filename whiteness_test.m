% This code is derived from the following sources:
% http://dsp.stackexchange.com/questions/7678/determining-the-whiteness-of-noise
% http://www.mathworks.com/matlabcentral/fileexchange/43312-tests-against-the-null-hyptothesis-of-whiteness
function [h, p_value, R] = whiteness_test (x, maxlag, alpha)
    x = bsxfun(@minus, x, mean(x, 2));

    N = size(x, 2);

    if nargin < 2
        maxlag = N - 1;
    end
    
    if nargin < 3
        alpha = 0.05;
    end
    
    num_channels = size(x, 1);
    h = zeros(num_channels, 1);
    p_value = zeros(num_channels, 1);
    R = zeros(num_channels, 1);
    
    for i=1:num_channels
        % 'biased'   - scales the raw cross-correlation by 1/M.
        % Estimates of covariance function
        [r, lag] = xcorr(x(i,:), maxlag, 'biased');
        
        % Test statistic
        R(i) = N ./ r(:, lag == 0).^2 .* sum(r(:, lag > 0).^2, 2);
        
        % If x is zero-mean white, then R is ~X^2(m).
        p_value(i) = 1 - chi2cdf(R(i), maxlag);
        T = chi2inv (1-alpha, maxlag);
        
        % h == 1 means non-white
        h(i) = R(i) > T;
    end
end