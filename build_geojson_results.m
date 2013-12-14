function xhats = build_geojson_results ()
    datasets = {
        'airportloop', ...
        'cudtrt13triphammercu', ...
        'dtrt13triphammer', ...
        'ncayugast', ...
        'rt13warrenrd', ...
    };

    static_gain_modes = {
        'static-slow',
        'static-fast',
    };

    kalman_gain_modes = {
        'kalman-steady-state', ...
        'kalman-time-varying', ...
    };

    process_gains = {0.1, 1, 10};
    adaptive_process_gains = {false, true};
    
    xhats = {};

    for i=1:length(datasets)
        dataset = datasets{i};
        disp([dataset])
        for j=1:length(static_gain_modes)
            gain_mode = static_gain_modes{j};
            disp(['=> ' gain_mode])
            xhat = filter_dataset(dataset, 5, 'feedback', gain_mode, 1, false);
            name = [dataset '-' strrep(gain_mode, '-', '_')];
            xhats{end+1} = [];
            xhats{end}.name = name;
            xhats{end}.xhat = xhat;
            save_as_geojson(['receiver_paths/' name '.json'], xhat, name, false);
        end
        for j=1:length(kalman_gain_modes)
            gain_mode = kalman_gain_modes{j};
            disp(['=> ' gain_mode])
            for k=1:length(process_gains)
                for l=1:length(adaptive_process_gains)
                    adaptive_process_gain = adaptive_process_gains{l};
                    process_gain = process_gains{k};
                    xhat = filter_dataset(dataset, 5, 'feedback', gain_mode, process_gain, adaptive_process_gain);
                    name = [dataset '-' strrep(gain_mode, '-', '_') '_' strrep(num2str(process_gain),'.','_') '_' num2str(adaptive_process_gain)];
                    xhats{end+1} = [];
                    xhats{end}.name = name;
                    xhats{end}.xhat = xhat;
                    save_as_geojson(['receiver_paths/' name '.json'], xhat, name, true);
                end
            end
        end
        xhat = filter_dataset(dataset, 5, 'absolute', 1, false);
        name = [dataset '-absolute'];
        xhats{end+1} = [];
        xhats{end}.name = name;
        xhats{end}.xhat = xhat;
        save_as_geojson(['receiver_paths/' name '.json'], xhat, name, false);
        xhat = filter_dataset(dataset, 5, 'integrate', 1, false);
        name = [dataset '-integrate'];
        xhats{end+1} = [];
        xhats{end}.name = name;
        xhats{end}.xhat = xhat;
        save_as_geojson(['receiver_paths/' name '.json'], xhat, name, false);
    end
end
