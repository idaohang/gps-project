function ephem = load_ephem(dataset_name)
    ephem = load(['ephem_' dataset_name '.asc']);
    
    if size(ephem,2) == 22
       ephem = [ephem,ephem(:,4)*[0 1]];
    end
end