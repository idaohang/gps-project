% findsat.m	(actual file name: findsat.m)
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
%	< INCLUDES CORRECTIONS >
%
% this GPS function computes satellite positions given the orbital
% ephemerides for at given times
%
% input: 'ephem' matrix which rows contain orbital ephemerides for
%		a given satellite
%					< see formatData.m for description >
%	 		't' vector conatins the GPS times each satellite position
%		will be calculated at; if all satellite positions are being
%		calculated at the same GPS time, 't' can be a scalar with
%		that GPS time; if only one satellite is specified in 'ephem'
%		and 't' contains a vector of times, then the satellite
%		position is found at multiple times
%
% output: 'satLoc' matrix which rows contain an SV id number, a GPS
%		time (seconds), and the ECEF coordinates (meters) of the
%		satellite location
%						[ svID GPStime ECEFx ECEFy ECEFz ;
%		  				  svID GPStime ECEFx ECEFy ECEFz ;
%											...
%						  svID GPStime ECEFx EFECy ECEFz ]
%
%	< INCLUDES CORRECTIONS >
%
function satLoc = findsat(ephem,t)
% define physical constants
constant;
% determine number of satellites; exit if zero satellites
satellites = size(ephem,1);
if (satellites == 0)
    return;
end
% define orbital parameters
t0e = ephem(:,4);		% ephemeris reference time (seconds)
ecc = ephem(:,5);		% eccentricity (unitless)
sqrta = ephem(:,6);		% square root of semi-major axis (meters1/2)
omega0 = ephem(:,7);	% argument of perigee (radians)
M0 = ephem(:,8);		% mean anomaly at reference time (radians)
l0 = ephem(:,9);		% right ascension at reference (radians)
lDot = ephem(:,10);     % rate of right acension (radians/second)
dn = ephem(:,11);		% mean motion difference (radians/second)
i0 = ephem(:,12); 		% inclination angle at reference time (radians)
iDot = ephem(:,13); 	% inclination angle rate (radians/second)
cuc = ephem(:,14);		% latitude cosine harmonic correction (radians)
cus = ephem(:,15);		% latitude sine harmonic correction (radians)
crc = ephem(:,16);		% orbit radius cosine harmonic correction (meters)
crs = ephem(:,17);		% orbit radius sine harmonic correction (meters)
cic = ephem(:,18);		% inclination cosine harmonic correction (radians)
cis = ephem(:,19);		% inclination sine harmonic correction (radians)
% if time parameter 't' is a vector and only one satellite; find
% that satellite position over an array of times
if ((size(t,1) ~= 1) && (satellites == 1))
    sat_ephem = ephem;
    for samples = 2:size(t,1)
        ephem = [ ephem ; sat_ephem ];
    end
end
% if time parameter 't' is a single value, create time vector
if (size(t,1) == 1)
    t = t .* ones(satellites,1);
end
% define time of position request and the change in t from epoch; 
% correct for possible week crossovers; 604800 seconds in a GPS week.
t_corr = t - 604800*round((t - t0e)*(1/604800));
dt = t_corr - t0e;
% calculate mean anomaly with corrections
n0 = sqrt(muearth) * (sqrta).^(-3);
n = n0 + dn;
M = M0 + n .* dt;
% compute the eccentric anomaly from mean anomaly using
% Newton-Raphson method to solve for 'E' in:
%		f(E) = M - E + ecc * sin(E) = 0
E = M;
for k = 1:10
    f = M - E + ecc .* sin(E);
    dfdE = - 1 + ecc .* cos(E);
    dE = - f ./ dfdE;
    E = E + dE;
end
% calculate true anomoly from eccentric anomoly
sinnu = sqrt(1 - ecc.^2) .* sin(E) ./ (1 - ecc .* cos(E));
cosnu = (cos(E) - ecc) ./ (1 - ecc .* cos(E));
nu = atan2(sinnu,cosnu);
% calculate the argument of latitude and use it to calculate the
% harmonic correction to the argument of perigee along with
% the total argument of perigee
Phi = nu + omega0;
twoPhi = 2*Phi;
cos2Phi = cos(twoPhi);
sin2Phi = sin(twoPhi);
omegaCorr = cuc .* cos2Phi + cus .* sin2Phi;
omega = omega0 + omegaCorr;
% calculate longitude of ascending node with secular correction
lcorr = lDot.*dt;
l = l0 - OmegaE .* t_corr + lcorr;
% calculate orbital radius with harmonic corrections
rCorr = crc .* cos2Phi + crs .* sin2Phi;
r = (sqrta.^2) .* (1 - ecc .* cos(E)) + rCorr;
% calculate inclination with secular and harminic corrections
iCorr = cic .* cos2Phi + cis .* sin2Phi + iDot .* dt;
i = i0 + iCorr;
% find position in orbital plane that has the x axis pointed
% along the ascending node and the z axis pointed along the
% orbital angular momentum.
u = omega + nu;
xp = r .* cos(u);
yp = r .* sin(u);
% find satellite position in ECEF coordinates,
% after having multiplied the vector [xp;yp;0]
% by the matrix Rz(-l)*Rx(-i) "by hand" in order to
% derive appropriate formulas for each of the 3 ECEF coordinates.
ECEFx = (xp .* cos(l)) - (yp .* cos(i) .* sin(l));
ECEFy = (xp .* sin(l)) + (yp .* cos(i) .* cos(l));
ECEFz = (yp .* sin(i));
% return satellite locations
satLoc = [ ephem(:,1) t ECEFx ECEFy ECEFz ];
return;