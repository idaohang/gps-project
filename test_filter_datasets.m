function test_filter_datasets()
    datasets = {
        'cudtrt13triphammercu',
        'dtrt13triphammer',
        'ncayugast',
        'rt13warrenrd',
    };

    gain_modes = {
        'static-slow',
        'static-fast',
        'kalman-steady-state',
        'kalman-time-varying',
    };

    for i=1:length(datasets)
        dataset = datasets{i};
        printf([dataset '\n'])
        for j=1:length(gain_modes)
            gain_mode = gain_modes{j};
            printf(['\t' gain_mode '\n'])
            xhat = filter_dataset(dataset, 'feedback', gain_mode);
            draw_trajectory(xhat);
            drawnow
            eigenvalues = plot_max_abs_eigenvalues(xhat);
            assert (max (eigenvalues) < 1)
            drawnow
            jerk = compute_jerk (xhat);
        end
        xhat = filter_dataset(dataset, 'absolute');
        draw_trajectory (xhat);
        drawnow
        xhat = filter_dataset(dataset, 'integrate');
        draw_trajectory (xhat);
        drawnow
    end
end