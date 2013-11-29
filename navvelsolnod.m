% navvelsolnod.m		(actual file name: navvelsolnod.m)
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% < OVERDETERMINED SATELLITE POSITION SOLUTION AND VELOCITY 
%   SOLUTION WITH IONOSPHERIC CORRECTIONS AND NEUTRAL 
%   ATMOSPHERE CORRECTIONS IN THE POSITION SOLUTION
%   AND WITH ENFORCEMENT OF ELEVATION MASK LOWER LIMIT >
%
% This program determines a navigational solution directly from raw
% pseudo-ranges and satellite ephemerides.  It also determines
% the velocity solution directly from raw carrier Doppler
% shifts.  The user is required to supply the name of the
% file that contains the Doppler shift data, an initial receiver 
% location (within 100 kilometers), four or more satellite SV ids  
% to use in the navigational solution from a list of available 
% SV ids, the pressure, temperature, and relative humidity at 
% the receiver, and a pseudo-range/Doppler-shift measurement 
% sample which specifies a GPS time, pseudo-range measurements,
% and carrier Doppler shift measurements to be used to calculate 
% the navigational solution.
%
% The files required for navvelsolnod.m are:
%
%		* constant.m
%       * deltl1ionocalc.m
%       * deltnacalc.m
%		* ecef.m
%       * elevazim.m
%		* findsatvelclock.m
%		* formatdata.m
%		* latlong.m
%		* printGPS.m
%		* solveposvelod.m
%		
% Ephemeris data, pseudo-range data, carrier Doppler shift data
% and ionosphere model data are needed for calculating the 
% corrected pseudo-ranges, the satellite positions, the
% satellite velocities, the final navigational solution,
% and the final velocity solution.  The program parseRinex.m 
% must have already formated data into four corresponding files; 
% ephem.asc, ion.asc, obs.asc (containing the L1 pseudorange
% data), and XXXXXXX.mat (containing the L1 carrier Doppler
% shift data in a Matlab workspace variable named obs). 
%
% navvelsolnod.m consists of the following steps and function calls:
%
%		* opens data files - ephem.asc, ion.asc, and obs.asc
%       * inputs name of beat carrier Doppler shift data file and
%         opens it and stores it in an appropriate array.
%       * inputs neutral atmosphere pressure (p), temparature (TdegC),
%         and relative humidity (hrel) data and computes temperature
%         in degrees Kelvin (TdegK).  If empty arrays are input, 
%         then nominal values will be used: p = 1013.25 millibars,
%         TdegK =  288.15 deg K, and hrel = 0.50.
%       * inputs flags that determine whether ionosphere delay
%         corrections and neutral atmosphere delay corrections
%         are applied during solution for position.
%		* inputs the observation station location guess
%       * inputs elevation mask angle.
%		* inputs four or more satellite SVs to be used
%		* calls formatdata to format 'ephem' and 'obs' data
%         structures
%		* inputs which pseudorange/Doppler-shift measurement is to
%         be used to calculate the position
%		* determine pseudorange measurements, carrier Doppler
%         shift measurements, and GPS time of sample chosen
%		* calls solveposvelod which iteratively calculates the receiver's
%         position from satellite locations and pseudoranges.
%         afterwards, it calculates the receiver's velocity
%         and its clock error rate.
%       * prints out satellite usage, elevation, and azimuth information
%         along with DOP and pseudorange RMS error information
%		* calls printGPS to print out results 
%       * inputs the true observation station location and
%         uses it to compute errors in vertical/east/north 
%         coordinates and display them
%       * Computes vertical/east/north velocity components
%         and displays them along with receiver clock error rate.
%
% < OVERDETERMINED SATELLITE POSITION SOLUTION AND VELOCITY 
%   SOLUTION WITH IONOSPHERIC CORRECTIONS AND NEUTRAL 
%   ATMOSPHERE CORRECTIONS IN THE POSITION SOLUTION
%   AND WITH ENFORCEMENT OF ELEVATION MASK LOWER LIMIT >
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
% retrieve pseudorange data from input file -- pseudo.asc
	load obs.asc; 
	pseudoData = obs;
	clear obs;
% retrieve carrier Doppler shift data from input file named by user.
    doppdatafilename = ...
         input('\nEnter name of carrier Doppler shift .mat file : ','s');
    load(doppdatafilename);
    DoppData = obs;
    clear obs
% check whether first sample, time, and SV columns of Doppdata and 
% pseudoData agree.  Quit with an error if they do not.
    if (size(pseudoData,1) ~= size(DoppData,1)) | ...
            (size(pseudoData,2) ~= size(DoppData,2))
       disp('Error in navvelsolnod.m: Doppler data array and')
       disp(' pseudorange data array have different dimensions.')
       return;
    end
    idum = [1,2,(3:2:size(pseudoData,2))]';
    if norm(pseudoData(:,idum) - DoppData(:,idum)) > 0
       disp('Error in navvelsolnod.m: Doppler data array and')
       disp(' pseudorange data array fail to have completely')
       disp(' identical sample indices, sample times,')
       disp(' and SV ID numbers.')
       return;
    end
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
       disp('Error in navvelsolnod.m: The following SVs')
       disp(' appear to lack valid pseudorange data:')
       fprintf('%d ',SVbadvec);
       disp(' ')
       return;
    end
	gpsTime = pseudo(s,1);
% retrieve the carrier Doppler shift data for the sample.
    nsatsuseddum = size(ephem,1);
    Doppshift = zeros(nsatsuseddum,1);
    SVs_Doppshift_raw = DoppData(s,3:2:end)';
    Doppshift_raw = DoppData(s,4:2:end)';
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
    clear nsatsuseddum SVs_Doppshift_raw Doppshift_raw k SVk idumk
    if size(SVbadvec,1) > 0
       disp(' ')
       disp('Error in navvelsolnod.m: The following SVs')
       disp(' appear to lack valid carrier Doppler shift data:')
       fprintf('%d ',SVbadvec);
       disp(' ')
       return;
    end
    clear SVbadvec
% call solveposvelod passing pseudoranges 'pseudoR', carrier
% Doppler shift Doppshift, an initial positional guess 'guess', 
% a GPS time 'gpsTime', and ephemeris data 'ephem', and other 
% relevant parameters.
    [posvelOBS,DOP,el_az,SVsused,sigmaPR,sigmaDopp] = ...
                   solveposvelod(ephem,pseudoR,Doppshift,guess,...
                                 gpsTime,ionParam,iflagion,elevmask,...
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
    fprintf('%9.4f ',DOP(:,1));
    disp(' ')
    disp(['RMS pseudorange residual, sigmaPR = ',...
           num2str(sigmaPR),' m.'])
    disp(' ')
    disp(['    GDOPv     PDOPv     TDOPv     HDOPv',...
          '     VDOPv'])
    fprintf('%9.4f ',DOP(:,2));
    disp(' ')
    disp(['RMS carrier Doppler shift residual, sigmaDopp = ',...
           num2str(sigmaDopp),' Hz.'])
% output the navigational solution
	printGPS(posvelOBS(1,1:5));
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
    ECEFerror = posvelOBS(1,2:4)' - obsPos_true;
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
    clear ECEFerror dx dy dz
    disp(' ')
    disp('Navigation Solution Errors :')
    disp('****************************')
	fprintf('\nVertical error        = %9.4f meters',Verror);
	fprintf('\nEast error            = %9.4f meters',Eerror);
	fprintf('\nNorth error           = %9.4f meters',Nerror);
    fprintf('\nTotal Error magnitude = %9.4f meters \n',MAGerror);
% compute velocity in vertical/east/north coordinates, display
% it, and display the receiver clock rate error solution.
    ECEFxdot = posvelOBS(1,6);
    ECEFydot = posvelOBS(1,7);
    ECEFzdot = posvelOBS(1,8);
	Vvel = cos(latOBS)*(ECEFxdot*cos(longOBS) + ...
		     ECEFydot*sin(longOBS)) + ECEFzdot*sin(latOBS);
	Evel = ECEFydot*cos(longOBS) - ECEFxdot*sin(longOBS);
	Nvel = ECEFzdot*cos(latOBS) - sin(latOBS)*...
		     (ECEFxdot*cos(longOBS) + ECEFydot*sin(longOBS));
    clear latOBS longOBS ECEFxdot ECEFydot ECEFzdot
    disp(' ')
    disp('Velocity Solution :')
    disp('*******************')
	fprintf('\nVertical velocity          = %9.4f meters/sec',Vvel);
	fprintf('\nEast velocity              = %9.4f meters/sec',Evel);
	fprintf('\nNorth velocity             = %9.4f meters/sec',Nvel);
    fprintf('\nReceiver clock rate error  = %d seconds/sec \n',...
             posvelOBS(1,9));    