function jerk = compute_jerk (x_hat)
    jerk = zeros(3, length(x_hat)-1);
    for i=1:length(x_hat)-1
        jerk(:,i) = (x_hat{i+1}.acceleration - x_hat{i}.acceleration)/(x_hat{i+1}.time - x_hat{i}.time);
    end
    lla = latlong(x_hat{1}.position');
    phi = lla(1)*pi/180;
    lambda = lla(2)*pi/180;
    A_VEN_ECEF = ...
      [cos(phi),0,sin(phi);0,1,0;-sin(phi),0,cos(phi)]*...
      [cos(lambda),sin(lambda),0;...
        -sin(lambda),cos(lambda),0;0,0,1];
    jerk = A_VEN_ECEF * jerk;
end
