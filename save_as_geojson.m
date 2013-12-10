function save_as_geojson (filename, x_hat, name)
    f = fopen(filename, 'wb');

    fprintf(f, '{\n');
    fprintf(f, '\t"type":"Feature",\n');
    fprintf(f, '\t"geometry": {\n');
    fprintf(f, '\t\t"type": "LineString",\n');
    fprintf(f, '\t\t"coordinates": [');

    for i=1:length(x_hat)
        if i > 1
            fprintf(f, ', ');
        end
        position = x_hat{i}.position;
        lla = latlong(position');
        fprintf(f, '[%f, %f]', lla(2), lla(1));
    end

    fprintf(f, ']\n');
    fprintf(f, '\t},\n');
    fprintf(f, '\t"properties": {\n');
    fprintf(f, '\t\t"name": "%s"\n', name);
    fprintf(f, '\t}\n');
    fprintf(f, '}\n');

    fclose(f);
end
