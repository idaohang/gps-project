function x_hat = filter_dataset(dataset_name, gain_mode)
    constant;

    dataset = load_dataset(dataset_name);
    
    SV = visible_satellite_filter(dataset.pseudorange(1,:));
    
    [ephem1 pseudo] = formatdata(dataset.ephem,dataset.pseudorange(1,:),SV');
    Doppshift = format_doppler_shift(ephem1, dataset.doppler_shift(1,:));

        
    pseudoR = pseudo(3:2:end)';
    guess = ecef([42.44 -76.48 200]);
    gpsTime = dataset.pseudorange(1,2);
    
    iflagna = true;
    iflagion = true;
    elevmask = 5;
    % TODO replace solveposod with solveposvelod
    [posOBS,~,~,~,~] = ...
                   solveposvelod(ephem1,pseudoR,Doppshift,guess,gpsTime,...
                              dataset.ion_params,iflagion,elevmask,dataset.weather.p,...
                              dataset.weather.TdegK,dataset.weather.hrel,iflagna);
                          
    num_samples = size(dataset.pseudorange,1);
    
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
    
    for k=1:num_samples-1
        t_Rk = dataset.pseudorange(k, 2);
        t_Rkp1 = dataset.pseudorange(k+1, 2);
        
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

        pseudorange_kp1 = dataset.pseudorange(k+1, :);
        doppler_shift_kp1 = dataset.doppler_shift(k+1, :);
        
        assert (pseudorange_kp1(1) == k+1);
        assert (doppler_shift_kp1(1) == k+1);
        
        SV = visible_satellite_filter(pseudorange_kp1);

        [ephem_kp1 pseudo] = formatdata(dataset.ephem,pseudorange_kp1,SV');
        Doppshift = format_doppler_shift(ephem_kp1, doppler_shift_kp1);
    
        pseudoR = pseudo(3:2:end)';

        [posOBS,~,~,~,~,~,Q,Qv] = ...
               solveposvelod_DOP(ephem_kp1,pseudoR,Doppshift, x_hat.position(:,k)',t_Rkp1,...
                          dataset.ion_params,iflagion,elevmask,dataset.weather.p,...
                          dataset.weather.TdegK,dataset.weather.hrel,iflagna);
          
        y_kp1 = [];
        y_kp1.time = posOBS(1);
        y_kp1.position = posOBS(2:4)';
        y_kp1.velocity = posOBS(6:8)';
        y_kp1.clock_offset = c*posOBS(5);
        y_kp1.clock_rate_offset = c*posOBS(9);
        
        v_kp1 = [];
        v_kp1.position = y_kp1.position - x_kp1_bar.position;
        v_kp1.velocity = y_kp1.velocity - x_kp1_bar.velocity;
        v_kp1.clock_offset = y_kp1.clock_offset - x_kp1_bar.clock_offset;
        v_kp1.clock_rate_offset = y_kp1.clock_rate_offset - x_kp1_bar.clock_rate_offset;
        
        x_k_hat = [];
        x_k_hat.time = x_hat.time(k);
        x_k_hat.position = x_hat.position(:,k);
        x_k_hat.velocity = x_hat.velocity(:,k);
        x_k_hat.acceleration = x_hat.acceleration(:,k);
        x_k_hat.clock_offset = x_hat.clock_offset(k);
        x_k_hat.clock_rate_offset = x_hat.clock_rate_offset(k);
        x_kp1_hat = apply_feedback_correction(x_kp1_bar, v_kp1, delta_T_Rk_clock_rate_offset, Q, Qv, y_kp1, x_k_hat, gain_mode);
        
        x_hat.time(k+1) = x_kp1_hat.time;
        x_hat.position(:,k+1) = x_kp1_hat.position;
        x_hat.velocity(:,k+1) = x_kp1_hat.velocity;
        x_hat.acceleration(:,k+1) = x_kp1_hat.acceleration;
        x_hat.clock_offset(k+1) = x_kp1_hat.clock_offset;
        x_hat.clock_rate_offset(k+1) = x_kp1_hat.clock_rate_offset;
    end
    

end