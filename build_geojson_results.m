function build_geojson_results ()
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

    parfor i=1:length(datasets)
        dataset = datasets{i};
        disp([dataset])
        for j=1:length(static_gain_modes)
            gain_mode = static_gain_modes{j};
            disp(['=> ' gain_mode])
            xhat = filter_dataset(dataset, 'feedback', gain_mode, 1);
            save_as_geojson(['receiver_paths/' dataset '-' strrep(gain_mode, '-', '_') '.json'], xhat, [dataset ' ' gain_mode], false);
        end
        for j=1:length(kalman_gain_modes)
            gain_mode = kalman_gain_modes{j};
            disp(['=> ' gain_mode])
            for k=1:length(process_gains)
                process_gain = process_gains{k};
                xhat = filter_dataset(dataset, 'feedback', gain_mode, process_gain);
                save_as_geojson(['receiver_paths/' dataset '-' strrep(gain_mode, '-', '_') '_' strrep(num2str(process_gain),'.','_') '.json'], xhat, [dataset ' ' gain_mode ' ' num2str(process_gain)], true);
            end
        end
        xhat = filter_dataset(dataset, 'absolute');
        save_as_geojson(['receiver_paths/' dataset '-absolute.json'], xhat, [dataset '-absolute'], false);
        xhat = filter_dataset(dataset, 'integrate');
        save_as_geojson(['receiver_paths/' dataset '-integrate.json'], xhat, [dataset '-integrate'], false);
    end
end
