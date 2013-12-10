function [names, xhats] = build_geojson_results ()
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
    
    xhats = {};
    names = {};

    for i=1:length(datasets)
        dataset = datasets{i};
        disp([dataset])
        for j=1:length(static_gain_modes)
            gain_mode = static_gain_modes{j};
            disp(['=> ' gain_mode])
            xhat = filter_dataset(dataset, 'feedback', gain_mode, 1);
            name = [dataset '-' strrep(gain_mode, '-', '_')];
            names{end+1} = name;
            xhats{end+1} = xhat;
            save_as_geojson(['receiver_paths/' name '.json'], xhat, name, false);
        end
        for j=1:length(kalman_gain_modes)
            gain_mode = kalman_gain_modes{j};
            disp(['=> ' gain_mode])
            for k=1:length(process_gains)
                process_gain = process_gains{k};
                xhat = filter_dataset(dataset, 'feedback', gain_mode, process_gain);
                name = [dataset '-' strrep(gain_mode, '-', '_') '_' strrep(num2str(process_gain),'.','_')];
                names{end+1} = name;
                xhats{end+1} = xhat;
                save_as_geojson(['receiver_paths/' name '.json'], xhat, name, true);
            end
        end
        xhat = filter_dataset(dataset, 'absolute');
        name = [dataset '-absolute'];
        names{end+1} = name;
        xhats{end+1} = xhat;
        save_as_geojson(['receiver_paths/' name '.json'], xhat, name, false);
        xhat = filter_dataset(dataset, 'integrate');
        name = [dataset '-integrate'];
        names{end+1} = name;
        xhats{end+1} = xhat;
        save_as_geojson(['receiver_paths/' name '.json'], xhat, name, false);
    end
end
