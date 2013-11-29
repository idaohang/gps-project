% rangecalc.m	(actual file name: rangecalc.m)
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% this function calculates the ranges from the observation station
% to selected satellites based on point-to-point distance
%
% input: 'ephem' matrix which rows contain orbital ephemerides for
%		a given satellite
%					< see formatData.m for description >
%			'pseudo' matrix which rows contain pseudo-range samples
%		for a given time
%					< see formatData.m for description >
%					< in this case, psuedo = obs structure >
%	 		'obsLoc' vector which contains the ECEF coordinates
%		(meters) of a reference position which ranges will be
%		calculated from
%						[ ECEFx ECEFy ECEFz ]
%
% output: 'range' matrix which rows contain a GPS time (seconds),
% 		and then pairs of SV id numbers with corresponding calculated
% 		ranges (meters)
%						[ GPStime svID r svID r ... ;
%						  GPStime svID r svID r ... ;
%											...
%						  GPStime svID r svID r ... ]
%
function range = rangecalc(ephem,pseudo,obsLoc)
% define physical constants
constant;
% clear variable 'range'
range = [ ];
% determine time samples
GPStime = pseudo(:,1);
% detemine number of samples taken
samples = size(GPStime,1);
% determine number of satellites ranges being calculated
satellites = size(ephem,1);
% define observation station location in ECEF coordinates
obsX = obsLoc(1);
obsY = obsLoc(2);
obsZ = obsLoc(3);
% create 'range' by calculating ranges from observation location to
% each satellite for each time sample
for t = 1:samples
  % calculate initial range based on GPS time sample
  satXYZ = findsat(ephem,GPStime(t));
  satX = satXYZ(:,3);
  satY = satXYZ(:,4);
  satZ = satXYZ(:,5);
  r = ((satX - obsX).^2 + (satY - obsY).^2 + ...
    (satZ - obsZ).^2).^0.5;
  % use this range 'r' to recalculate satellite locations based
  % on transmission time = reception time - range/c, which is a
  % satellite movement correction
  transT = r ./ c;
  satXYZ = findsat(ephem,GPStime(t) - transT);
  satX = satXYZ(:,3);
  satY = satXYZ(:,4);
  satZ = satXYZ(:,5);
  % correction due to rotation of the Earth, also based 
  % based here in the ad hoc range/c delay calculation.
  deltheta = transT .* OmegaE;
  cosdeltheta = cos(deltheta);
  sindeltheta = sin(deltheta);
  satX_uncorr = satX;
  satY_uncorr = satY;
  satX =   satX_uncorr .* cosdeltheta + satY_uncorr .* sindeltheta;
  satY = - satX_uncorr .* sindeltheta + satY_uncorr .* cosdeltheta;
  % calculate actual satellite range from observation station
  r = ((satX - obsX).^2 + (satY - obsY).^2 + ...
    (satZ - obsZ).^2).^0.5;
  % add range calculations into 'range' structure
  sample = GPStime(t);
  for i = 1:satellites
    sample = [ sample ephem(i,1) r(i) ];
  end
  range = [ range; sample ];
end
% return calculated ranges
return