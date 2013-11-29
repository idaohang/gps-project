function SV = visible_satellite_filter (obs)    
    SV = obs(1,3:2:end);
    SV = unique(SV(:));
    if SV(1,1) == 0
       SV(1,:) = [];
    end
    nsatsdum = size(SV,1);
    for k = nsatsdum:-1:1
       if ~any(SV(k,1) == SV(:,1))
          SV(k,:) = [];
       end
    end
end