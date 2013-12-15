function dataset = load_dataset(dataset_name)
    dataset = [];

    original_dataset_name = dataset_name;
    if strcmp(dataset_name, 'airportloop_outage') || ...
       strcmp(dataset_name, 'airportloop_corruption')
        dataset_name = 'airportloop';
    end
    
    dataset.pseudorange = load(['obs_' dataset_name '.asc']);
    dataset.doppler_shift = load(['obsdopp_' dataset_name '.asc']);
    dataset.ion_params = load(['ion_' dataset_name '.asc']);
    dataset.ephem = load_ephem(dataset_name);
    dataset.weather = load_weather(dataset_name);
    
    if strcmp(original_dataset_name, 'airportloop_outage')
        dataset.pseudorange(100:105,3:end) = 0;
        dataset.doppler_shift(100:105,3:end) = 0;
    end

    if strcmp(original_dataset_name, 'airportloop_corruption')
        for i=100:105
            for j=3:2:size(dataset.pseudorange,2)
                disp(dataset.pseudorange(i,j))
                if dataset.pseudorange(i,j) == 11
                    dataset.pseudorange(i,j+1) = dataset.pseudorange(i,j+1) + 15^2*randn(1);
                end
            end
        end
    end
    
    assert (size(dataset.pseudorange,1) == size(dataset.doppler_shift,1));
end
