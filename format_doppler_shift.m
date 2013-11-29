function Doppshift = format_doppler_shift(ephem, obsdopp)
    nsatsuseddum = size(ephem,1);
    Doppshift = zeros(nsatsuseddum,1);
    SVs_Doppshift_raw = obsdopp(1,3:2:end)';
    Doppshift_raw = obsdopp(1,4:2:end)';
    SVbadvec = [];
    for k = 1:nsatsuseddum
       SVk = ephem(k,1);
       idumk = find(ephem(k,1) == SVs_Doppshift_raw);
       if size(idumk,1) < 1
          SVbadvec = [SVbadvec;SVk];
       else
          Doppshift(k,1) = Doppshift_raw(idumk(1,1),1);
       end
    end
    assert (size(SVbadvec, 1) == 0);
end