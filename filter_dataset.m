function x_hat = filter_dataset(dataset_name, mode, gain_mode, process_gain)
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
    [posOBS,~,~,~,sigmaPR,sigmaDopp] = ...
                   solveposvelod(ephem1,pseudoR,Doppshift,guess,gpsTime,...
                              dataset.ion_params,iflagion,elevmask,dataset.weather.p,...
                              dataset.weather.TdegK,dataset.weather.hrel,iflagna);
                          
    if isnan(sigmaPR)
        sigmaPR = 2;
    end
    
    if isnan(sigmaDopp)
        sigmaDopp = 0.15;
    end
                          
    num_samples = size(dataset.pseudorange,1);
    
    x_hat = cell(num_samples, 1);

    x_hat{1} = [];
    x_hat{1}.time = posOBS(1) - posOBS(5);
    x_hat{1}.position = posOBS(2:4)';
    x_hat{1}.velocity = posOBS(6:8)';
    x_hat{1}.acceleration = zeros(3, 1);
    x_hat{1}.clock_offset = c*posOBS(5);
    x_hat{1}.clock_rate_offset = c*posOBS(9);
    x_hat{1}.covariance = eye(11);
    x_hat{1}.Phi = zeros(11);
    
    for k=1:num_samples-1
        t_Rk = x_hat{k}.time;
        t_Rkp1 = dataset.pseudorange(k+1, 2);
        
        pseudorange_kp1 = dataset.pseudorange(k+1, :);
        doppler_shift_kp1 = dataset.doppler_shift(k+1, :);
        
        SV = visible_satellite_filter(pseudorange_kp1);

        [ephem_kp1 pseudo] = formatdata(dataset.ephem,pseudorange_kp1,SV');
        Doppshift = format_doppler_shift(ephem_kp1, doppler_shift_kp1);
    
        pseudoR = pseudo(3:2:end)';

        [posOBS,~,~,~,sigmaPRk,sigmaDoppk,Q,Qv] = ...
               solveposvelod_DOP(ephem_kp1,pseudoR,Doppshift, x_hat{k}.position',t_Rkp1,...
                          dataset.ion_params,iflagion,elevmask,dataset.weather.p,...
                          dataset.weather.TdegK,dataset.weather.hrel,iflagna);
                      
        if ~isnan(sigmaPRk)
            sigmaPR = sigmaPRk;
        end
        
        if ~isnan(sigmaDoppk)
            sigmaDopp = sigmaDoppk;
        end

        measurements = [];
        measurements.time = posOBS(1) - posOBS(5);
        measurements.position = posOBS(2:4)';
        measurements.velocity = posOBS(6:8)';
        measurements.clock_offset = c*posOBS(5);
        measurements.clock_rate_offset = c*posOBS(9);
        measurements.sigmaPR = sigmaPR;
        measurements.sigmaDopp = sigmaDopp;
        measurements.Q = Q;
        measurements.Qv = Qv;

        dt = measurements.time - t_Rk;
                
        deltaTRk = dt/(1 + x_hat{k}.clock_rate_offset/c);
       

        if strcmp(mode, 'feedback')
            x_hat{k+1} = feedback(deltaTRk, x_hat{k}, measurements, gain_mode, process_gain);
        elseif strcmp(mode, 'absolute')
            x_hat{k+1} = measurements;
        elseif strcmp(mode, 'integrate')
            x_hat{k+1} = measurements;
            x_hat{k+1}.position = x_hat{k}.position + measurements.velocity;
        end
    end
    

end