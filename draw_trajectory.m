function draw_trajectory(x_hat)
    mean_ecef = mean(x_hat.position, 2);
    mean_lla = latlong(mean_ecef');
    phiavg = mean_lla(1,1)*pi/180;
    lambdaavg = mean_lla(1,2)*pi/180;
    A_VEN_ECEF = ...
       [cos(phiavg),0,sin(phiavg);0,1,0;-sin(phiavg),0,cos(phiavg)]*...
       [cos(lambdaavg),sin(lambdaavg),0;...
          -sin(lambdaavg),cos(lambdaavg),0;0,0,1];
    local_position = A_VEN_ECEF*x_hat.position;
    velocity_mag = sqrt(sum(x_hat.velocity.^2));
    scatter(local_position(2,:), local_position(3,:), 10, velocity_mag);
    axis equal
end