function [mh, sh, mv, sv] = error_statistics (estimated, correct)
    assert(length(estimated) == length(correct));
    len = length(estimated);
    error = zeros(3, len);
    for i=1:len
        error(:, i) = estimated{i}.position - correct{i}.position;
    end
    lla = latlong(correct{1}.position');
    phi = lla(1)*pi/180;
    lambda = lla(2)*pi/180;
    A_VEN_ECEF = ...
      [cos(phi),0,sin(phi);0,1,0;-sin(phi),0,cos(phi)]*...
      [cos(lambda),sin(lambda),0;...
        -sin(lambda),cos(lambda),0;0,0,1];
    error = A_VEN_ECEF * error;
    mv = mean(error(1,:));
    sv = sqrt(var(error(1,:)));
    mh = mean(sqrt(sum(error(2:3,:).^2)));
    sh = sqrt(var(sqrt(sum(error(2:3,:).^2))));
    plot(1:size(error, 2), error(1,:));
    plot(1:size(error, 2), sqrt(sum(error(2:3,:).^2)));
end