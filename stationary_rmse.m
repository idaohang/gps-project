function [rmse, errors] = stationary_rmse (xhat, surveyed)
    positions = zeros(3, length(xhat));
    for i=1:length(xhat)
        positions(:,i) = xhat{i}.position;
    end

    diffs = bsxfun(@minus, positions, surveyed');

    errors = sqrt(sum(diffs.^2));

    rmse = sqrt(mean(errors.^2))
end
