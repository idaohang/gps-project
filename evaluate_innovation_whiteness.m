function retval = evaluate_innovation_whiteness(x_hat)
    innovations = zeros(8, length(x_hat)-1);
    lla = latlong(x_hat{1}.position');
    phi = lla(1)*pi/180;
    lambda = lla(2)*pi/180;
    A_VEN_ECEF = ...
      [cos(phi),0,sin(phi);0,1,0;-sin(phi),0,cos(phi)]*...
      [cos(lambda),sin(lambda),0;...
        -sin(lambda),cos(lambda),0;0,0,1];

    for i=2:length(x_hat)
        if false%isfield(x_hat{i}.metadata, 'S')
            innovations(:,i) = x_hat{i}.metadata.S^(-0.5) * x_hat{i}.innovation;
        else
            innovations(:,i) = x_hat{i}.innovation;
        end
    end
    innovations(1:3,:) = A_VEN_ECEF * innovations(1:3,:);
    innovations(4:6,:) = A_VEN_ECEF * innovations(4:6,:);

    h = whiteness_test_portmanteau(innovations(4:6,:), 20, 0.01)
end
