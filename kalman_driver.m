function kalman_driver(dataset_name)
    constant;

    obs = load(['obs_' dataset_name '.asc']);
    obsdopp = load(['obsdopp_' dataset_name '.asc']);
    ion = load(['ion_' dataset_name '.asc']);
    ephem = load_ephem(dataset_name);
    weather = load_weather(dataset_name);
    
    assert (size(obs,1) == size(obsdopp,1));
    
    SV = visible_satellite_filter(obs(1,:));
    
    [ephem1 pseudo] = formatdata(ephem,obs(1,:),SV');
    Doppshift = format_doppler_shift(ephem1, obsdopp(1,:));

        
    pseudoR = pseudo(3:2:end)';
    guess = ecef([42.44 -76.48 200]);
    gpsTime = obs(1,2);
    
    iflagna = true;
    iflagion = true;
    elevmask = 5;
    % TODO replace solveposod with solveposvelod
    [posOBS,~,~,~,~] = ...
                   solveposvelod(ephem1,pseudoR,Doppshift,guess,gpsTime,...
                              ion,iflagion,elevmask,weather.p,...
                              weather.TdegK,weather.hrel,iflagna);
                          
    num_samples = size(obs,1);
    
    x_hat = struct( ...
        'time', zeros(1, num_samples), ...
        'position', zeros(3, num_samples), ...
        'velocity', zeros(3, num_samples), ...
        'acceleration', zeros(3, num_samples), ...
        'clock_offset', zeros(1, num_samples), ...
        'clock_rate_offset', zeros(1, num_samples));

    x_hat.time(1) = posOBS(1);
    x_hat.position(:,1) = posOBS(2:4)';
    x_hat.velocity(:,1) = posOBS(6:8)';
    x_hat.acceleration(:,1) = zeros(3, 1);
    x_hat.clock_offset(1) = c*posOBS(5);
    x_hat.clock_rate_offset(1) = c*posOBS(9);
    
    for k=1:size(obs,1)-1
        t_Rk = obs(k, 2);
        t_Rkp1 = obs(k+1, 2);
        
        dt = t_Rkp1 - t_Rk;
                
        delta_T_Rk_clock_rate_offset = dt/(1 + x_hat.clock_rate_offset(k)/c);

        x_kp1_bar = [];
        x_kp1_bar.time = t_Rkp1;
        x_kp1_bar.position = x_hat.position(:,k) + ...
            delta_T_Rk_clock_rate_offset * x_hat.velocity(:,k) + ...
            0.5 * delta_T_Rk_clock_rate_offset^2 * x_hat.acceleration(:,k);
        x_kp1_bar.velocity = x_hat.velocity(:,k) + ...
            delta_T_Rk_clock_rate_offset * x_hat.acceleration(:,k);
        x_kp1_bar.acceleration = x_hat.acceleration(:,k);
        x_kp1_bar.clock_offset = x_hat.clock_offset(k) + ...
            delta_T_Rk_clock_rate_offset * x_hat.clock_rate_offset(k);
        x_kp1_bar.clock_rate_offset = x_hat.clock_rate_offset(k);

        obs_kp1 = obs(k+1, :);
        obsdopp_kp1 = obsdopp(k+1, :);
        
        assert (obs_kp1(1) == k+1);
        assert (obsdopp_kp1(1) == k+1);
        
        SV = visible_satellite_filter(obs_kp1);

        [ephem_kp1 pseudo] = formatdata(ephem,obs_kp1,SV');
        Doppshift = format_doppler_shift(ephem_kp1, obsdopp_kp1);
    
        pseudoR = pseudo(3:2:end)';

        [posOBS,~,~,~,~] = ...
               solveposvelod(ephem_kp1,pseudoR,Doppshift, x_hat.position(:,k)',t_Rkp1,...
                          ion,iflagion,elevmask,weather.p,...
                          weather.TdegK,weather.hrel,iflagna);
          
        y_kp1 = [];
        y_kp1.time = posOBS(1);
        y_kp1.position = posOBS(2:4)';
        y_kp1.velocity = posOBS(6:8)';
        y_kp1.clock_offset = c*posOBS(5);
        y_kp1.clock_rate_offset = c*posOBS(9);
        
        %x_kp1_bar.position = y_kp1.position;
        %x_kp1_bar.velocity = y_kp1.velocity;
        %x_kp1_bar.clock_offset = y_kp1.clock_offset;
        %x_kp1_bar.clock_rate_offset = y_kp1.clock_rate_offset;
        
        v_kp1 = [];
        v_kp1.position = y_kp1.position - x_kp1_bar.position;
        v_kp1.velocity = y_kp1.velocity - x_kp1_bar.velocity;
        v_kp1.clock_offset = y_kp1.clock_offset - x_kp1_bar.clock_offset;
        v_kp1.clock_rate_offset = y_kp1.clock_rate_offset - x_kp1_bar.clock_rate_offset;
        
        x_kp1_hat = apply_feedback_correction(x_kp1_bar, v_kp1, delta_T_Rk_clock_rate_offset);
        
        x_hat.time(k+1) = x_kp1_hat.time;
        x_hat.position(:,k+1) = x_kp1_hat.position;
        x_hat.velocity(:,k+1) = x_kp1_hat.velocity;
        x_hat.acceleration(:,k+1) = x_kp1_hat.acceleration;
        x_hat.clock_offset(k+1) = x_kp1_hat.clock_offset;
        x_hat.clock_rate_offset(k+1) = x_kp1_hat.clock_rate_offset;
    end
    
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