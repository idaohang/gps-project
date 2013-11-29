% sattraj.m (actual file name: sattraj.m)
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% This program plots the trajectories of a set of visible satellites
% for a given time range. The satellite positions are found in ECEF 
% coordinates. Also, elevation-azimuth look angles to each satellite
% are computed given the location of an observation station.
%
% The files required for sattraj.m are:     * constant.m
%                                           * ecef.m
%                                           * elevazim.m 
%                                           * findsat.m    
%                                           * formatdata.m
%                                           * gmttogps.m
%                                           * latlong.m
%                                           * plotsat.m
%
% Ephemeris data is needed for calculating the satellite positions.
% The program gps.exe will format the data into the file 
% 'ephem.asc'.
%
% The sattraj.m consists of the following steps and function calls:
%       * opens data file - ephem.asc
%       * inputs GPS time to start satellite trajectories and inputs
%           the length of the time for the trajectories
%       * inputs satellite SVs to calculate trajectories for
%       * input observation station location
%       * calls formatData to format 'ephem' data structure
%       * loops findsat which finds satellite locations in ECEF
%           coordinates for a given time
%       * calls elevazim which converts ECEF coordinates into
%           elevation and azimuth
%       * calls plotSat to plot trajectories
%
%
% retrieve ephemeris data from input file -- ephem.asc
    load ephem.asc;
    ephemData = ephem; 
    clear ephem;
% input GMT time to compute position at and convert to GPS time
    fprintf('\nEnter date to begin trajectory plot at \n');
    date = input('  in the form of   "mmddyy"  :  ','s');
    fprintf('Enter GMT time to begin trajectory plot at \n');
    time = input('  in the form of   "hhmmss"  :  ','s');
    t = gmttogps(date,time);
   fprintf('\nGPS Time for trajectory start  :  %11.3f seconds \n', t);
% input the duration of the trajectory plot
    fprintf('\nEnter the duration for the trajectory plot \n');
    duration = input('  in units of minutes :  ');
    sRate = 3;          % sample rate
    samples = fix(duration / sRate);    
% determine which satellites a position will be found for
    fprintf('\nEnter the satellites to calculate positions for \n');
    fprintf('\nfrom the list of SV ids --> ');
   fprintf('%d ',ephemData(:,1));
   fprintf('\n');
    SV = input('    in the form of   "[ sv sv . . . ]"  :  ');
    % default is all satellites
    if (isempty(SV))
        SV = linspace(1,32,32);
    end
% input the observation station location in latitude-longitude
% coordinates 
    fprintf('\nEnter the location of the observation station \n');
    obsLoc = input('    in the form "[ lat long alt ]"  :  ');
    obsLoc = ecef(obsLoc);
% call formatData which reformats 'ephemData' into 'ephem'
    ephem = formatdata(ephemData,[ ],SV);
% find locations of satellites along trajectory
    elevation_azimuth = [];
    for s = 1:samples
        % find satellite ECEF coordinates for a given time
        satLoc = findsat(ephem,t + (s - 1) * sRate * 60);
      % find elevation and azimuth of the satellites
      el_az = elevazim(satLoc,obsLoc);
        % add current elevations and azimuth to past results
        elevation_azimuth = [elevation_azimuth; el_az];
    end
% plot the computed satellite trajectories
    plotsat(elevation_azimuth);
