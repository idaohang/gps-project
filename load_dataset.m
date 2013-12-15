function dataset = load_dataset(dataset_name)
    dataset = [];

    dataset.pseudorange = load(['obs_' dataset_name '.asc']);
    dataset.doppler_shift = load(['obsdopp_' dataset_name '.asc']);
    dataset.ion_params = load(['ion_' dataset_name '.asc']);
    dataset.ephem = load_ephem(dataset_name);
    dataset.weather = load_weather(dataset_name);
    
    assert (size(dataset.pseudorange,1) == size(dataset.doppler_shift,1));

    dataset.pseudorange(1,2)
end
