function build_geojson_results ()
    datasets = {
        'cudtrt13triphammercu', ...
        'dtrt13triphammer', ...
        'ncayugast', ...
        'rt13warrenrd', ...
    };

    gain_modes = {
        'static-slow', ...
        'static-fast', ...
        'kalman-steady-state', ...
        'kalman-time-varying', ...
    };

    parfor i=1:length(datasets)
        dataset = datasets{i};
        disp([dataset])
        for j=1:length(gain_modes)
            gain_mode = gain_modes{j};
            disp(['=> ' gain_mode])
            xhat = filter_dataset(dataset, 'feedback', gain_mode);
            save_as_geojson(['receiver_paths/' dataset '-' strrep(gain_mode, '-', '_') '.json'], xhat, [dataset '-' gain_mode]);
        end
        xhat = filter_dataset(dataset, 'absolute');
        save_as_geojson(['receiver_paths/' dataset '-absolute.json'], xhat, [dataset '-absolute']);
        xhat = filter_dataset(dataset, 'integrate');
        save_as_geojson(['receiver_paths/' dataset '-integrate.json'], xhat, [dataset '-integrate']);
    end
end