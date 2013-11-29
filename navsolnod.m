% navsolnod.m		(actual file name: navsolnod.m)
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% < OVERDETERMINED SATELLITE SOLUTION WITH IONOSPHERIC CORRECTIONS
%   AND NEUTRAL ATMOSPHERE CORRECTIONS AND WITH ENFORCEMENT OF
%   ELEVATION MASK LOWER LIMIT >
%
% This program determines a navigational solution directly from raw
% pseudo-ranges and satellite ephemerides. The user is required to 
% supply an initial receiver location (within 100 kilometers), four
% or more satellite SV ids to use in the navigational solution 
% from a list of available SV ids, the pressure, temperature, 
% and relative humidity at the receiver, and a pseudo-range 
% measurement sample which specifies a GPS time and pseudo-range 
% measurements to be used to calculate the navigational solution.
%
% The files required for navsolnod.m are:
%
%		* constant.m
%       * deltl1ionocalc.m
%       * deltnacalc.m
%		* ecef.m
%       * elevazim.m
%		* findsatclock.m
%		* formatdata.m
%		* latlong.m
%		* printGPS.m
%		* solveposod.m
%		
% Ephemeris data, pseudo-range data, and ionosphere model data
% are needed for calculating the corrected pseudo-ranges, the 
% satellite positions, and the final navigational solution. 
% The program parseRinex.m must have already formated data into  
% three corresponding files; ephem.asc, ion.asc, and obs.asc. 
%
% navsolnod.m consists of the following steps and function calls:
%
%		* opens data files - ephem.asc, ion.asc, and obs.asc
%       * inputs neutral atmosphere pressure (p), temparature (TdegC),
%         and relative humidity (hrel) data and computes temperature
%         in degrees Kelvin (TdegK).  If empty arrays are input, 
%         then nominal values will be used: p = 1013.25 millibars,
%         TdegK =  288.15 deg K, and hrel = 0.50.
%       * inputs flags that determine whether ionosphere delay
%         corrections and neutral atmosphere delay corrections
%         are applied during solution of navigation problem.
%		* inputs the observation station location guess
%       * inputs elevation mask angle.
%		* inputs four or more satellite SVs to be used
%		* calls formatdata to format 'ephem' and 'obs' data
%         structures
%		* inputs which pseudo-range measurement is to be used
%         to calculate the position
%		* determine pseudo-range measurements and GPS time of
%         sample chosen
%		* calls solveposod which iteratively calculates the receiver's
%         position from satellite locations and pseudo-ranges 
%       * prints out satellite usage, elevation, and azimuth information
%         along with DOP and pseudorange RMS error information
%		* calls printGPS to print out results 
%       * inputs the true observation station location and
%         uses it to compute errors in vertical/east/north 
%         coordinates.
%
% < OVERDETERMINED SATELLITE SOLUTION WITH IONOSPHERIC CORRECTIONS
%   AND NEUTRAL ATMOSPHERE CORRECTIONS AND WITH ENFORCEMENT OF
%   ELEVATION MASK LOWER LIMIT >
%

% clear the Matlab workspace
    clear
% retrieve ephemeris data from input file -- ephem.asc 
	load ephem.asc;
	ephemData = ephem; 
	clear ephem;
% retrieve ionosphere model data from input file -- ion.asc
    load ion.asc
    ionParam = ion;
    clear ion
% retrieve pseudo-range data from input file -- pseudo.asc
	load obs.asc; 
	pseudoData = obs;
	clear obs; 
% input neutral atmosphere data.  Note: these data may be
% retrieved from a web site such as http://www.wunderground.com.
% this web site includes history data for past days.  Note that
% there are 33.86 millbars per inch of Mercury (in Hg)
	fprintf('\nEnter the atmospheric pressure at the receiver in \n');    
	p = input('millibars (also hPa), p : ');
	fprintf('\nEnter the temperature at the receiver in \n');    
	TdegC = input('degrees Centigrade, TdegC : ');
    TdegK = TdegC + 273.15;
    clear TdegC;
	fprintf('\nEnter the relative humidity at the receiver as \n');    
	hrel = input('a fraction in the range 0 to 1, hrel : ');
% input the decisions of whether to use the ionosphere delay
% and the neutral atmosphere delay in calculating the receiver
% clock error and the pseudorange errors.
	fprintf('\nEnter a "1" if ionospheric delay corrections are to \n');    
	iflagion = input('be used.  Enter a "0" otherwise : ');
	fprintf('\nEnter a "1" if neutral atmosphere delay corrections \n');    
	iflagna = input('are to be used.  Enter a "0" otherwise : ');    
% input initial guess of receiver's position in latitude-longitude
% coordinates 
	fprintf('\nEnter the approximate location of the observation station \n');
	guess = input(' in the form "[ latitude longitude altitude ]"  :  ');
    guess = ecef(guess);
% input elevation mask angle. 
	fprintf('\nEnter the elevation mask angle below which ');
	fprintf('\n the correspoonding pseudorange will not \n');
	elevmask = input(' be used (deg) elevmask :  ');
% get the list of SVIDs for which obs and ephem data are available.
    SVIDlist = pseudoData(:,3:2:end);
    SVIDlist = unique(SVIDlist(:));
    if SVIDlist(1,1) == 0
       SVIDlist(1,:) = [];
    end
    nsatsdum = size(SVIDlist,1);
    for k = nsatsdum:-1:1
       if ~any(SVIDlist(k,1) == ephemData(:,1))
          SVIDlist(k,:) = [];
       end
    end
    clear nsatsdum;
% determine which satellites will be used in the navigational
% solution
	fprintf('\nEnter the satellites to be used in the navigational solution ');
	fprintf('\nEnter FOUR or more - they must have pseudoranges:');
	fprintf('\nfrom the list of SV ids --> ');
    fprintf('%d ',SVIDlist);
    fprintf('\n');
	SV = input('	in the form of   "[ sv sv . . . ]"  :  ');
% call formatData which will reformat 'ephemData' and 'pseudoData'
% into the structures 'ephem' and 'pseudo' respectively
	[ ephem pseudo ] = formatdata(ephemData,pseudoData,SV);
% check that four satellites have been specified
	if (size(ephem,1) < 4)
		fprintf(['\nNeed to specify four or more satellites',...
                 ' with ephemerides !!!\n\n']);
		return;
	end
% determine which pseudo-range measurement will be used
	fprintf('\nEnter the pseudo-range measurement sample, which range');
	fprintf('\nfrom 1 to %d, to be used to calculate the position \n',...
            size(pseudo,1));
	s = input('	in the form of  "sample #"    :  ');
% check that valid sample number has been entered
	if ((s < 1) || (s > size(pseudo,1)))
		fprintf('\nSample number is out of range !!!\n\n');
		return;
	end
% determine pseudo-range measurements 'pseudoR' and GPS time 
% 'gpsTime' of the sample chosen
	pseudoR = pseudo(s,3:2:end)';
% check that have psuedo-range data for all satellites for the 
% sample
    if any(pseudoR == 0)
       idumbadvec = find(pseudoR == 0);
       SVbadvec = pseudo(s,2:2:end)';
       SVbadvec = SVbadvec(idumbadvec,1);
       disp(' ')
       disp('Error in navsolnod.m: The following SVs')
       disp(' appear to lack valid pseudorange data:')
       fprintf('%d ',SVbadvec);
       disp(' ')
       return;
    end
	gpsTime = pseudo(s,1);
% call solveposed passing pseudo-ranges 'pseudoR', an initial 
% positional guess 'guess', a GPS time 'gpsTime', and ephemeris 
% data 'ephem', and other relevant parameters.
    [posOBS,DOP,el_az,SVsused,sigmaPR] = ...
                   solveposod(ephem,pseudoR,guess,gpsTime,...
                              ionParam,iflagion,elevmask,...
                              p,TdegK,hrel,iflagna);
% output information about the solution.
    Nsatsused = size(SVsused,1);
    el_azused = zeros(Nsatsused,4);
    for jj = 1:Nsatsused
       idumjj = find(SVsused(jj,1) == el_az(:,1));
       el_azused(jj,:) = el_az(idumjj,:);
    end
    clear jj idumjj
    disp(' ')
    disp('Used satellites, elevations, & azimuths :')
    disp('*****************************************')
    disp(' ')
    disp('SV    elev (deg)   azim (deg)')
    fprintf('%2i    %7.2f      %7.2f\n',el_azused(:,[1 3 4])');
    disp(' ')
    disp('Solution accuracy parameters :')
    disp('******************************')
    disp(' ')
    disp(['    GDOP      PDOP      TDOP      HDOP',...
          '      VDOP'])
    fprintf('%9.4f ',DOP);
    disp(' ')
    disp(['RMS pseudorange residual, sigmaPR = ',...
           num2str(sigmaPR),' m.'])
% output the navigational solution
	printGPS(posOBS);
% input receiver's true position in latitude-longitude
% coordinates
    if iflagion == 0
       disp('Ionosphere delay corrections not used.')
    end
    if iflagna == 0
       disp('Neutral atmosphere delay corrections not used.')
    end
	fprintf('\nEnter the true location of the observation station \n');
	obsPoslatlongalt_true = ...
        input(['   in the form "[ latitude(deg) longitude(deg)',...
               ' altitude(m) ]"  :  ']);
    obsPos_true = ecef(obsPoslatlongalt_true);
    obsPos_true = obsPos_true(:);
% compute the error in ECEF coordinates and then transform
% it to vertical/east/north coordinates.
    ECEFerror = posOBS(1,2:4)' - obsPos_true;
    constant
    latOBS = obsPoslatlongalt_true(1)*degrad;		% radians
    longOBS = obsPoslatlongalt_true(2)*degrad;		% radians
    clear AA OmegaE degrad f0 lambdaL1 muearth BB c esquare...
          fL1 leapSeconds
    dx = ECEFerror(1,1);
    dy = ECEFerror(2,1);
    dz = ECEFerror(3,1);
	Verror = cos(latOBS)*(dx*cos(longOBS) + ...
		     dy*sin(longOBS)) + dz*sin(latOBS);
	Eerror = dy*cos(longOBS) - dx*sin(longOBS);
	Nerror = dz*cos(latOBS) - sin(latOBS)*...
		     (dx*cos(longOBS) + dy*sin(longOBS));
    MAGerror = norm(ECEFerror);
    clear ECEFerror latOBS longOBS dx dy dz
    disp(' ')
    disp('Navigation Solution Errors :')
    disp('****************************')
	fprintf('\nVertical error        = %9.4f meters',Verror);
	fprintf('\nEast error            = %9.4f meters',Eerror);
	fprintf('\nNorth error           = %9.4f meters',Nerror);
    fprintf('\nTotal Error magnitude = %9.4f meters \n',MAGerror);