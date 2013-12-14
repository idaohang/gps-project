function draw_trajectory(xhat)
    positions = cell2mat(cellfun(@(x) x.position, xhat, 'UniformOutput', false)');
    velocities = cell2mat(cellfun(@(x) x.velocity, xhat, 'UniformOutput', false)');
    mean_ecef = mean(positions, 2);
    mean_lla = latlong(mean_ecef');
    phiavg = mean_lla(1,1)*pi/180;
    lambdaavg = mean_lla(1,2)*pi/180;
    A_VEN_ECEF = ...
       [cos(phiavg),0,sin(phiavg);0,1,0;-sin(phiavg),0,cos(phiavg)]*...
       [cos(lambdaavg),sin(lambdaavg),0;...
          -sin(lambdaavg),cos(lambdaavg),0;0,0,1];
    positions = bsxfun(@minus, positions, mean_ecef);
    local_position = A_VEN_ECEF*positions;
    velocity_mag = sqrt(sum(velocities.^2));
    scatter(local_position(2,3:end), local_position(3,3:end), 10, velocity_mag(3:end));
    axis equal
end