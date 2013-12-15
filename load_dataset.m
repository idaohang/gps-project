function dataset = load_dataset(dataset_name)
    dataset = [];

    original_dataset_name = dataset_name;
    if strcmp(dataset_name, 'airportloop_outage')
        dataset_name = 'airportloop';
    end
    
    dataset.pseudorange = load(['obs_' dataset_name '.asc']);
    dataset.doppler_shift = load(['obsdopp_' dataset_name '.asc']);
    dataset.ion_params = load(['ion_' dataset_name '.asc']);
    dataset.ephem = load_ephem(dataset_name);
    dataset.weather = load_weather(dataset_name);
    
    if strcmp(original_dataset_name, 'airportloop_outage')
        %dataset.pseudorange(100:100,3:end) = 0;
        %dataset.doppler_shift(100:100,3:end) = 0;
        dataset.pseudorange(100:105,3:end) = 0;
        dataset.doppler_shift(100:105,3:end) = 0;
        %dataset.pseudorange = dataset.pseudorange(374:395,:);
        %dataset.doppler_shift = dataset.doppler_shift(374:395,:);
    end
    
    assert (size(dataset.pseudorange,1) == size(dataset.doppler_shift,1));
end
