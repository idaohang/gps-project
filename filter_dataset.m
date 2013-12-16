function x_hat = filter_dataset(dataset_name, elevmask, bias, mode, gain_mode, process_gain, measurement_gain)
    constant;

    dataset = load_dataset(dataset_name);
    
    SV = visible_satellite_filter(dataset.pseudorange(1,:));
    
    [ephem1, pseudo] = formatdata(dataset.ephem,dataset.pseudorange(1,:),SV');
    Doppshift = format_doppler_shift(ephem1, dataset.doppler_shift(1,:));
        
    pseudoR = pseudo(3:2:end)';
    guess = ecef([42.44 -76.48 200]);
    gpsTime = dataset.pseudorange(1,2);
    
    iflagna = true;
    iflagion = true;
    [posOBS,~,~,~,sigmaPR,sigmaDopp] = ...
                   solveposvelod(ephem1,pseudoR,Doppshift,guess,gpsTime,...
                              dataset.ion_params,iflagion,elevmask,dataset.weather.p,...
                              dataset.weather.TdegK,dataset.weather.hrel,iflagna);
                          
    num_samples = size(dataset.pseudorange,1);
    
    pseudorange_bias_history = zeros(32, num_samples);
    doppler_shift_bias_history = zeros(32, num_samples);
    sigma_pr_history = zeros(1, num_samples);
    sigma_dopp_history = zeros(1, num_samples);
    t_history = zeros(1, num_samples);
    
    if isnan(sigmaPR)
        sigma_pr_history(1) = 2;
    else
        sigma_pr_history(1) = sigmaPR;
    end
    
    if isnan(sigmaDopp)
        sigma_dopp_history(1) = 0.15;
    else
        sigma_dopp_history(1) = sigmaDopp;
    end
                          
    x_hat = cell(num_samples, 1);

    if bias
        num_states = 75;
    else
        num_states = 11;
    end
    
    x_hat{1} = [];
    x_hat{1}.time = posOBS(1) - posOBS(5);
    x_hat{1}.position = posOBS(2:4)';
    x_hat{1}.velocity = posOBS(6:8)';
    x_hat{1}.acceleration = zeros(3, 1);
    x_hat{1}.clock_offset = c*posOBS(5);
    x_hat{1}.clock_rate_offset = c*posOBS(9);
    x_hat{1}.covariance = 1e6*eye(num_states);
    x_hat{1}.Phi = zeros(num_states);
    x_hat{1}.innovation = zeros(8, 1);
    x_hat{1}.process_gain = process_gain;
    x_hat{1}.gain = zeros(11,8);
    x_hat{1}.gain(1:6,1:6) = eye(6);
    x_hat{1}.gain(10:11,7:8) = eye(2);
    x_hat{1}.t = 0;
    if bias
        x_hat{1}.pseudorange_bias = zeros(32, 1);
        x_hat{1}.doppler_shift_bias = zeros(32, 1); 
    end
    
    window = 1;    
    for k=1:num_samples-1
        t_Rk = x_hat{k}.time;
        t_Rkp1 = dataset.pseudorange(k+1, 2);
        
        pseudorange_kp1 = dataset.pseudorange(k+1, :);
        doppler_shift_kp1 = dataset.doppler_shift(k+1, :);
        
        SV = visible_satellite_filter(pseudorange_kp1);
        
        if isempty(SV)
            A = process_matrix (deltaTRk, bias);
            x = state_vector (x_hat{k});
            Ax = A*x;
            Q = zeros(11);
            Q(7:9,7:9) = process_gain*eye(3);
            P = A*x_hat{k}.covariance*A' + Q;
            x_hat{k+1} = create_state (Ax, t_Rkp1, zeros(8,1), P, [], [], [], [], 0);
            continue
        end

        [ephem_kp1, pseudo] = formatdata(dataset.ephem,pseudorange_kp1,SV');
        Doppshift = format_doppler_shift(ephem_kp1, doppler_shift_kp1);
    
        pseudoR = pseudo(3:2:end)';
        pseudorange_bias = zeros(size(pseudoR));
        doppler_shift_bias = zeros(size(Doppshift));
        if bias
            for i=1:length(SV)
                if SV(i) <= 32
                    pseudorange_bias(i) = x_hat{k}.pseudorange_bias(SV(i));
                    doppler_shift_bias(i) = x_hat{k}.doppler_shift_bias(SV(i));
                end
            end
        end

        [posOBS,~,~,SVsused,sigmaPRk,sigmaDoppk,Q,Qv,pseudorange_bias, doppler_shift_bias] = ...
               solveposvelod_DOP(ephem_kp1,pseudoR,Doppshift, x_hat{k}.position',t_Rkp1,...
                          dataset.ion_params,iflagion,elevmask,dataset.weather.p,...
                          dataset.weather.TdegK,dataset.weather.hrel,iflagna, pseudorange_bias, doppler_shift_bias);
        if bias 
            [~,~,~,~,~,~,~,~,pseudorange_bias, doppler_shift_bias] = ...
                   solveposvelod_DOP(ephem_kp1,pseudoR,Doppshift, x_hat{k}.position',t_Rkp1,...
                              dataset.ion_params,iflagion,elevmask,dataset.weather.p,...
                              dataset.weather.TdegK,dataset.weather.hrel,iflagna);
        end
          
        if ~isnan(sigmaPRk)
            sigma_pr_history(k+1) = sigmaPRk;
        end
        
        if ~isnan(sigmaDoppk)
            sigma_dopp_history(k+1) = sigmaDoppk;
        end

        measurements = [];
        measurements.time = posOBS(1) - posOBS(5);
        measurements.position = posOBS(2:4)';
        measurements.velocity = posOBS(6:8)';
        measurements.clock_offset = c*posOBS(5);
        measurements.clock_rate_offset = c*posOBS(9);

        sigma_window = 60;
        sigma_pr_data = sigma_pr_history(max(1,k-sigma_window):k+1);
        sigma_dopp_data = sigma_dopp_history(max(1,k-sigma_window):k+1);
        measurements.sigmaPR = sqrt(median(sigma_pr_data(sigma_pr_data ~= 0).^2));
        measurements.sigmaDopp = sqrt(median(sigma_dopp_data(sigma_dopp_data ~= 0).^2));
       
        measurements.Q = Q;
        measurements.Qv = Qv;
        
        bias_window = 20;
        if bias
            if ~isempty(pseudorange_bias)
                pseudorange_var = zeros(size(pseudorange_bias));
                doppler_shift_var = zeros(size(doppler_shift_bias));
                for i=1:length(SVsused)
                    pseudorange_bias_history(SVsused(i),k+1) = pseudorange_bias(i);
                    doppler_shift_bias_history(SVsused(i),k+1) = doppler_shift_bias(i);
                
                    if k > bias_window
                        pr_data = pseudorange_bias_history(SVsused(i),k-bias_window:k+1);
                        pr_data = pr_data(pr_data ~= 0);
                        dop_data = doppler_shift_bias_history(SVsused(i),k-bias_window:k+1);
                        dop_data = dop_data(dop_data ~= 0);
                        pseudorange_bias(i) = pr_data(end);
                        doppler_shift_bias(i) = dop_data(end);
                        pseudorange_var(i) = 1e5*var(pr_data(pr_data ~= 0));
                        doppler_shift_var(i) = 1e5*var(dop_data(dop_data ~= 0));
                    else
                        pseudorange_bias(i) = 0;
                        doppler_shift_bias(i) = 0;
                        pseudorange_var(i) = 0;
                        doppler_shift_var(i) = 0;
                    end
                end
                measurements.pseudorange_bias = [SVsused pseudorange_bias];
                measurements.doppler_shift_bias = [SVsused doppler_shift_bias];
                measurements.pseudorange_var = [SVsused pseudorange_var];
                measurements.doppler_shift_var = [SVsused doppler_shift_var];
            else
                measurements.pseudorange_bias = zeros(0,2);
                measurements.doppler_shift_bias = zeros(0,2);
                measurements.pseudorange_var = zeros(0,2);
                measurements.doppler_shift_var = zeros(0,2);
            end
        end

        dt = measurements.time - t_Rk;
                
        deltaTRk = dt/(1 + x_hat{k}.clock_rate_offset/c);
        
        if strcmp(mode, 'feedback')
            x_hat{k+1} = feedback(deltaTRk, x_hat{k}, measurements, gain_mode, process_gain, measurement_gain);
            if isfield(x_hat{k+1}.metadata, 'S')
                S = x_hat{k+1}.metadata.S;
                v = x_hat{k+1}.innovation;
                t = v(1:8)' / S(1:8,1:8) * v(1:8);
                t_history(k+1) = t;
                t = t_history(t_history ~= 0);
                if length(t) >= 20
                    t = sum(t)/length(t);
        
                    n_v = 8;
                    N = length(t);

                    alpha = 0.05;
                    r1 = chi2inv( alpha/2,  N*n_v ) / N;
                    r2 = chi2inv( 1 - alpha/2,  N*n_v ) / N;
                    pass = (t > r1) && (t < r2);
                    pass = true;
                    if ~pass
                        disp('fail')
                        A = process_matrix (deltaTRk, bias);
                        x = state_vector (x_hat{k});
                        Ax = A*x;
                        Q = zeros(11);
                        Q(7:9,7:9) = process_gain*eye(3);
                        P = A*x_hat{k}.covariance*A' + Q;
                        x_hat{k+1} = create_state (Ax, measurements.time, zeros(8,1), P, [], [], [], [], 0);
                        t_history(k+1) = 0;
                    end
                end
            end
        elseif strcmp(mode, 'absolute')
            x_hat{k+1} = measurements;
        elseif strcmp(mode, 'integrate')
            x_hat{k+1} = measurements;
            x_hat{k+1}.position = x_hat{k}.position + deltaTRk*measurements.velocity;
        end
    end
end
