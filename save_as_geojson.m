function save_as_geojson (filename, x_hat, name, save_covariance)
    f = fopen(filename, 'wb');

    fprintf(f, '{\n');
    fprintf(f, '\t"type":"Feature",\n');
    fprintf(f, '\t"geometry": {\n');
    fprintf(f, '\t\t"type": "MultiLineString",\n');
    fprintf(f, '\t\t"coordinates": [\n');
    
    fprintf(f, '\t\t\t[');

    for i=1:length(x_hat)
        if i > 1
            fprintf(f, ', ');
        end
        position = x_hat{i}.position;
        lla = latlong(position');
        fprintf(f, '[%f, %f]', lla(2), lla(1));
    end

    fprintf(f, ']');
        
    if save_covariance
        num_covariances = min([length(x_hat), 1000]);
        sample_spacing = length(x_hat) / num_covariances;
        m = chi2inv(0.99, 3);
        
        
        for i=1:num_covariances
            covariance = x_hat{floor(sample_spacing * i)}.covariance;
            
            center = x_hat{floor(sample_spacing * i)}.position;
            lla = latlong(center');
            phi = lla(1)*pi/180;
            lambda = lla(2)*pi/180;
            A_VEN_ECEF = ...
              [cos(phi),0,sin(phi);0,1,0;-sin(phi),0,cos(phi)]*...
              [cos(lambda),sin(lambda),0;...
                -sin(lambda),cos(lambda),0;0,0,1];

           
            covariance_local = A_VEN_ECEF * covariance(1:3,1:3) * A_VEN_ECEF';
            
            covariance_EN = covariance_local(2:3, 2:3);
            
            fprintf(f, ',\n');
            fprintf(f, '\t\t\t[');
            
            [eigvec,eigval] = eig(m*covariance_EN);
            p=0:pi/16:2*pi;
            xy = [cos(p'),sin(p')] * sqrt(eigval) * eigvec';
            offset_global = A_VEN_ECEF' * [zeros(size(p)); xy'];
            
            for j=1:length(offset_global)
                position = center + offset_global(:,j);
                if j > 1
                    fprintf(f, ', ');
                end
                lla = latlong(position');
                fprintf(f, '[%f, %f]', lla(2), lla(1));
            end
            fprintf(f, ']');
        end
        fprintf(f, '\n');
    end
    fprintf(f, '\t\t]\n');
    
    fprintf(f, '\t},\n');
    fprintf(f, '\t"properties": {\n');
    fprintf(f, '\t\t"name": "%s"\n', name);
    fprintf(f, '\t}\n');
    fprintf(f, '}\n');

    fclose(f);
end
