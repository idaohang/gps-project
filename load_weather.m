function weather = load_weather (dataset_name)
    if strcmp(dataset_name, 'airportloop') || ...
       strcmp(dataset_name, 'ncayugast') || ...
       strcmp(dataset_name, 'rt13warrenrd')
       pinHg = 30.2;
       p = 33.86*pinHg; % convert to millibars (also hPa)
       TdegF = 53;
       TdegC = (TdegF - 32)*(5/9);
       TdegK = TdegC + 273.15;
       hrel = 0.94;
    elseif strcmp(dataset_name, 'dtrt13triphammer') || ...
           strcmp(dataset_name, 'cudtrt13triphammercu')
       pinHg = 30.3;
       p = 33.86*pinHg; % convert to millibars (also hPa)
       TdegF = 40;
       TdegC = (TdegF - 32)*(5/9);
       TdegK = TdegC + 273.15;
       hrel = 0.81; 
    elseif strcmp(dataset_name, 'stationary')
       p = 1015.46;
       TdegK = 4.83 + 273.15;
       hrel = 0.76;
    end
    
    weather = [];
    weather.p = p;
    weather.TdegK = TdegK;
    weather.hrel = hrel;
end