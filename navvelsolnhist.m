% navvelsolnhist.m		(actual file name: navvelsolnhist.m)
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% < FULL TIME HISTORY OVERDETERMINED SATELLITE POSITION SOLUTION  
%   AND VELOCITY SOLUTION WITH IONOSPHERIC CORRECTIONS AND NEUTRAL 
%   ATMOSPHERE CORRECTIONS IN THE POSITION SOLUTION,
%   WITH ENFORCEMENT OF ELEVATION MASK LOWER LIMIT, AND WITH
%   VELOCITY-SMOOTHED POSIITION OUTPUTS IN VERTICAL/EAST/
%   NORTH COORDINATES.>
%
% This program determines a navigational solution directly from raw
% pseudo-ranges and satellite ephemerides.  It also determines
% the velocity solution directly from raw carrier Doppler
% shifts.  It does this for an entire time history of data that
% is user-selectable from one of 3 supplied sets of GPS data files
% or from user-supplied files.  In the latter case, the user 
% is required to supply the name of the file that contains the 
% Doppler shift data.  The user also supplies an initial receiver 
% location (within 100 kilometers) and four or more satellite SV ids  
% to use in the navigational solution from a list of available 
% SV ids.  If the user provides his or her own input data,
% then he or she must also provide the pressure, temperature, 
% and relative humidity at the receiver.  In the case
% of the 3 supplied sets of data, the meteorological data
% are provided by the present script.
%
% This script also performs analysis of the results in
% vertical/east/north components.  It compares velocity-smoothed
% position (i.e., trapezoidally integrated velocity 
% with the bias removed beteen its average and the pseudorange
% solution time history's average) with pseudorange position
% in vertical/east/north coordinates.
%
% The files required for navvelsolnhist.m are:
%
%		* constant.m
%       * deltl1ionocalc.m
%       * deltnacalc.m
%		* ecef.m
%       * elevazim.m
%		* findsatvelclock.m
%		* formatdata.m
%       * integrate_vel_ven.m
%		* latlong.m
%		* printGPS.m
%		* solveposvelod.m
%
%       * ephem_ncayugast.asc
%       * obs_ncayugast.asc
%       * obsdopp_ncayugast.asc
%       * ion_ncayugast.asc
%
%       * ephem_rt13warrenrd.asc
%       * obs_rt13warrenrd.asc
%       * obsdopp_rt13warrenrd.asc
%       * ion_rt13warrenrd.asc
%
%       * ephem_airportloop.asc
%       * obs_airportloop.asc
%       * obsdopp_airportloop.asc
%       * ion_airportloop.asc
%
%       * ephem_dtrt13triphammer.asc
%       * obs_dtrt13triphammer.asc
%       * obsdopp_dtrt13triphammer.asc
%       * ion_dtrt13triphammer.asc
%
%       * ephem_cudtrt13triphammercu.asc
%       * obs_cudtrt13triphammercu.asc
%       * obsdopp_cudtrt13triphammercu.asc
%       * ion_cudtrt13triphammercu.asc
%
% The first 3 sets of files were collected by Mark Psiaki and
% Brent Ledvina using a Magellan receiver on 17 Oct. 2000 while
% driving a car around the Ithaca area.
%
% The data in the files
% obs_ncayugast.asc and obsdopp_ncayugast.asc constitute
% only the first 120 samples of the the data in the raw
% files obs_ncayugast_all.asc and obsdopp_ncayugast_all.asc.
% The last 53 samples have been deleted due to problems with
% sample spacings not being even at about 1 sec or due to
% increased pseudorange RMS residuals.
%
% The last two sets of files were collected by SK Callanan 
% on 1 Dec. 2011.  They were sent to Mark Psiaki on 31 Oct. 
% 2012.  The raw pseudorange and
% Doppler shift observables were sent in the file 
% obsprcpdopp_skdrive_01dec2011.asc. which has an unusual
% record format of gpsTime, SV_a, pseudorange_a, beatcarrierphase_c,
% Dopplershift_a, SV_b, pseudorange_b, beatcarrierphase_b,
% Dopplershift_b, ...
%
% These data were split into standard
% obs.asc-type files, obs_skdrive_01dec2011_all.asc for
% the pseudoranges and obsdopp_skdrive_01dec2011_all.asc for
% the Doppler shift data.  Their record formats are
% sampleNumber, gpsTime, SV_a, pseudorange_a, SV_b, pseudorange_b, ...
% for the pseudorange file and sampleNumber, gpsTime, SV_a, 
% Dopplershift_a, SV_b, Dopplershift_b, ... for the Doppler
% shift file.  The final observables files for the first data 
% set, obs_dtrt13triphammer.asc and obsdopp_dtrt13triphammer.asc 
% cover samples 187 through 1396 of the original data, but with a 
% total of 10 samples deleted this range of data, some because 
% they had fewer than 4 satellites with ephemerides in 
% ephem_skdrive_01dec2011.asc, which equals 
% ephem_dtrt13triphammer.asc, and some because they had
% a pseudorange RMS residual greater
% than 30 m, a GDOP greater than 30, or both.  This data set
% looks like a drive starting downtown going east and then west
% on a downtown street, then driving up Rt. 13 to Triphammer road, 
% followed by a drive down Triphammer Rd. about to the Cornell campus.
% 
% The final observables files for the of the second SK Callanan 
% data set, obs_cudtrt13triphammercu.asc and 
% obsdopp_cudtrt13triphammercu.asc cover samples 1397 
% through 3148 of the original data, but with a total 
% of 17 samples deleted this range of data, some because 
% they had fewer than 4 satellites with ephemerides in 
% ephem_skdrive_01dec2011.asc, which equals 
% ephem_cudtrt13triphammercu.asc, and some because they had
% a pseudorange RMS residual greater
% than 30 m, a GDOP greater than 30, or both.  This data set
% looks like a drive starting on or near the Cornell campus,
% proceding downtown, looping from downtown up Rt. 13
% to Triphammer Rd., and coming back to the Cornell campus
% on Triphammer Rd.
%		
% Ephemeris data, pseudo-range data, carrier Doppler shift data
% and ionosphere model data are needed for calculating the 
% corrected pseudo-ranges, the satellite positions, the
% satellite velocities, the final navigational solution,
% and the final velocity solution.  In the case of
% provided data, the needed information is in the selected 
% set of .asc files listed ablve.  In the case of user-
% supplied data, the program parseRinex.m must have already 
% been used to formate data into four corresponding files; 
% ephem.asc, ion.asc, obs.asc (containing the L1 pseudorange
% data), and XXXXXXX.mat (containing the L1 carrier Doppler
% shift data in a Matlab workspace variable named obs). 
%
% navvelsolnhist.m consists of the following steps and function calls:
%
%       * inputs decision of user about which data files to use
%		* opens data files, prompting user for name of Doppler
%         shift file if user-supplied data files have been selected.
%       * initializes or inputs neutral atmosphere pressure (p), 
%         temparature (TdegC), and relative humidity (hrel) data 
%         and computes temperature in degrees Kelvin (TdegK).  If 
%         empty arrays are input for a user-supplied data case, 
%         then nominal values will be used: p = 1013.25 millibars,
%         TdegK =  288.15 deg K, and hrel = 0.50.
%       * sets flags that specify that the ionosphere delay
%         corrections and neutral atmosphere delay corrections
%         will be applied during solutions of navigation solution.
%		* inputs the initial receiver location guess.
%       * inputs elevation mask angle.
%       * creates a receiver sample time history vector.
%       * sets up pseudorange and carrier Doppler shift data
%         arrays and array that keeps track of which satellites
%         are available at each sample time.  
%		* inputs four or more satellite SVs to be used
%         and prunes away data for additional SVs that
%         have not been selected.
%       * initializes arrays for storing time history results.
%       * cycles through all data points, solves for position
%         and velocity using solveposvelod.m, and stores results.
%       * plots north position vs. east position relative to average
%         point, both pure pseudorange-based solutions and
%         solution from time integral of velocity with
%         unknown bias removed by forcing mean solution
%         to match that ov pseudorange-based solution
%       * plots vertical position time history, both
%         from raw pseudoranges and from integrated 
%         velocity soltion with unknown bias removed, 
%         and plots clock error solution time history, 
%         both from raw ppseudoranges and from integrated 
%         receiver clock error rate soltion with Doppler 
%         removed.
%       * plots velocity solution time history and approximate
%         velocity solution created by finite-differencing position
%         solutions.
%       * plots GDOP, HDOP, and VDOP time histories.
%       * plots pseudorange and Doppler shift RMS error
%         time histories.
%
% < OVERDETERMINED SATELLITE POSITION SOLUTION AND VELOCITY 
%   SOLUTION WITH IONOSPHERIC CORRECTIONS AND NEUTRAL 
%   ATMOSPHERE CORRECTIONS IN THE POSITION SOLUTION
%   AND WITH ENFORCEMENT OF ELEVATION MASK LOWER LIMIT >
%

% clear the Matlab workspace
    clear
% input flag that tells which set of data files to use.
    disp(' ')
    fprintf('\nChoose the data set to use by entering:');    
    fprintf('\n   1 for North Cayuga St. run, 17 Oct. 2000');    
    fprintf('\n   2 for Rt. 13/Warren Rd. run, 17 Oct. 2000');    
    fprintf('\n   3 for Tompkins Airport Loops run, 17 Oct. 2000');    
    fprintf('\n   4 for Downtown/Rt. 13/Triphammer run, 01 Dec. 2011');    
    fprintf('\n   5 for CU/Downtown/Rt. 13/Triphammer/CU run, 01 Dec. 2011');    
    fprintf('\n   6 for user supplied files (i.e., ephem.asc,');    
    iflagfiles = ...
      input('\n     obs.asc, ion.asc, and Doppler .mat-file) : ');
% retrieve data files, ephemerides, pseudorange, carrier Doppler
% shift, and ionosphere model parameters.
    if iflagfiles == 1
       load ephem_ncayugast.asc
       ephemData = ephem_ncayugast;
       clear ephem_ncayugast
       load obs_ncayugast.asc
       pseudoData = obs_ncayugast;
       clear obs_ncayugast
       load obsdopp_ncayugast.asc
       DoppData = obsdopp_ncayugast;
       clear obsdopp_ncayugast
       load ion_ncayugast.asc
       ionParam = ion_ncayugast;
       clear ion_ncayugast
    elseif iflagfiles == 2
       load ephem_rt13warrenrd.asc
       ephemData = ephem_rt13warrenrd;
       clear ephem_rt13warrenrd
       load obs_rt13warrenrd.asc
       pseudoData = obs_rt13warrenrd;
       clear obs_rt13warrenrd
       load obsdopp_rt13warrenrd.asc
       DoppData = obsdopp_rt13warrenrd;
       clear obsdopp_rt13warrenrd
       load ion_rt13warrenrd.asc
       ionParam = ion_rt13warrenrd;
       clear ion_rt13warrenrd
    elseif iflagfiles == 3
       load ephem_airportloop.asc
       ephemData = ephem_airportloop;
       clear ephem_airportloop
       load obs_airportloop.asc
       pseudoData = obs_airportloop;
       clear obs_airportloop
       load obsdopp_airportloop.asc
       DoppData = obsdopp_airportloop;
       clear obsdopp_airportloop
       load ion_airportloop.asc
       ionParam = ion_airportloop;
       clear ion_airportloop
    elseif iflagfiles == 4
       load ephem_dtrt13triphammer.asc
       ephemData = ephem_dtrt13triphammer;
       clear ephem_dtrt13triphammer
       load obs_dtrt13triphammer.asc
       pseudoData = obs_dtrt13triphammer;
       clear obs_dtrt13triphammer
       load obsdopp_dtrt13triphammer.asc
       DoppData = obsdopp_dtrt13triphammer;
       clear obsdopp_dtrt13triphammer
       load ion_dtrt13triphammer.asc
       ionParam = ion_dtrt13triphammer;
       clear ion_dtrt13triphammer
    elseif iflagfiles == 5
       load ephem_cudtrt13triphammercu.asc
       ephemData = ephem_cudtrt13triphammercu;
       clear ephem_cudtrt13triphammercu
       load obs_cudtrt13triphammercu.asc
       pseudoData = obs_cudtrt13triphammercu;
       clear obs_cudtrt13triphammercu
       load obsdopp_cudtrt13triphammercu.asc
       DoppData = obsdopp_cudtrt13triphammercu;
       clear obsdopp_cudtrt13triphammercu
       load ion_cudtrt13triphammercu.asc
       ionParam = ion_cudtrt13triphammercu;
       clear ion_cudtrt13triphammercu
    else
       load ephem.asc;
       ephemData = ephem; 
       clear ephem;
       load obs.asc
       pseudoData = obs;
       clear obs
       % retrieve carrier Doppler shift data from input file named by user.
       doppdatafilename = ...
         input('\nEnter name of carrier Doppler shift .mat file : ','s');
       load(doppdatafilename);
       DoppData = obs;
       clear obs
       load ion.asc
       ionParam = ion;
       clear ion
    end
% add missing ephemData data if there are only 22 columns in ephemData.
% this expansion assumes that the clock correction reference time is
% the same as the ephemeris reference time and that the single
% frequency clock calibration bias value is zero.
    if size(ephemData,2) == 22
       ephemData = [ephemData,ephemData(:,4)*[0 1]];
    end
% check whether first sample, time, and SV columns of Doppdata and 
% pseudoData agree.  Quit with an error if they do not.
    if (size(pseudoData,1) ~= size(DoppData,1)) | ...
            (size(pseudoData,2) ~= size(DoppData,2))
       disp('Error in navvelsolnhist.m: Doppler data array and')
       disp(' pseudorange data array have different dimensions.')
       return;
    end
    idum = [1,2,(3:2:size(pseudoData,2))]';
    if norm(pseudoData(:,idum) - DoppData(:,idum)) > 0
       disp('Error in navvelsolnhist.m: Doppler data array and')
       disp(' pseudorange data array fail to have completely')
       disp(' identical sample indices, sample times,')
       disp(' and SV ID numbers.')
       return;
    end
% initialize or input neutral atmosphere data.  these data were
% retrieved from a web site such as http://www.wunderground.com.
% for 17 Oct. 2000 at 9 a.m. in Ithaca, NY: pressure of 30.2
% inches of Mercury, temperature of 53 degF, and relative
% humidity of 94% (choose average value for that day).
    if iflagfiles <= 3
% These data were
% retrieved from a web site such as http://www.wunderground.com.
% for 17 Oct. 2000 at 9 a.m. in Ithaca, NY: pressure of 30.2
% inches of Mercury, temperature of 53 degF, and relative
% humidity of 94% (choose average value for that day).
       pinHg = 30.2;
       p = 33.86*pinHg; % convert to millibars (also hPa)
       TdegF = 53;
       TdegC = (TdegF - 32)*(5/9);
       TdegK = TdegC + 273.15;
       hrel = 0.94;
       clear pinHg TdegF TdegC
    elseif (iflagfiles == 4) | (iflagfiles == 5)
% These data were
% retrieved from a web site such as http://www.wunderground.com.
% for 01 Dec. 2011 at 11 a.m. in Ithaca, NY: pressure of 30.3
% inches of Mercury, temperature of 40 degF, and relative
% humidity of 81% (choose average value for that day).
       pinHg = 30.3;
       p = 33.86*pinHg; % convert to millibars (also hPa)
       TdegF = 40;
       TdegC = (TdegF - 32)*(5/9);
       TdegK = TdegC + 273.15;
       hrel = 0.81;
       clear pinHg TdegF TdegC
    else
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
    end
% set the decisions to use the ionosphere delay
% and the neutral atmosphere delay in calculating the receiver
% clock error and the pseudorange errors.
	iflagion = 1;
	iflagna = 1;    
% input initial guess of receiver's position in latitude-longitude
% coordinates 
	fprintf('\nEnter the approximate initial location of the receiver \n');
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
    for jj = nsatsdum:-1:1
       if ~any(SVIDlist(jj,1) == ephemData(:,1))
          SVIDlist(jj,:) = [];
       end
    end
    clear jj nsatsdum
% determine number of samples and the receiver clock GPS times
    N = size(pseudoData,1);
    gpsTimehist = pseudoData(:,2);
% create a pure pseudorange data array, a Doppler shift data
% array, and an indexing array for the
% available pseudorange data and Doppler shift data.  
% Pseudorange and Doppler shift
% data in the arrays pseudoDatahist and
% DopplerDatahist are valid only if the corresponding
% entries of iPRNflagshist contain 1 values.  Each
% column corresponds to the GPS satellite whose SVID is
% contained in the corresponding row of SVIDlist, and each
% row corresponds to the time in the corresponding row
% of gpsTimehist.
    Nsats = size(SVIDlist,1);
    iPRNflagshist = zeros(N,Nsats);
    pseudoDatahist = zeros(N,Nsats);
    DopplerDatahist = zeros(N,Nsats);
    for k = 1:N
       SVIDlistk = pseudoData(k,3:2:end)';
       for jj = 1:Nsats
          idumkjj = find(SVIDlistk == SVIDlist(jj,1));
          if size(idumkjj,1) == 1
             iPRNflagshist(k,jj) = 1;
             idumkjj_mod = 2*(idumkjj + 1);
             pseudoDatahist(k,jj) = pseudoData(k,idumkjj_mod);
             DopplerDatahist(k,jj) = DoppData(k,idumkjj_mod);
          end
       end
    end
    clear k jj idumkjj idumkjj_mod pseudoData DoppData SVIDlistk
% determine number of times each SV is unavailable
    Notavailablelist = zeros(Nsats,1);
    for jj = 1:Nsats
       Notavailablelist(jj,1) = sum(iPRNflagshist(:,jj) == 0);
    end
% determine which satellites will be used in the navigational
% solution and in the velocity solution
	fprintf('\nEnter the satellites to be used in the navigational solution ');
	fprintf('\n Enter FOUR or more - they must have pseudoranges -- from ');
	fprintf('\n the list of SV ids (& numbers of samples unavailable out of %i) \n',N);
    fprintf('%4i ',SVIDlist);
    fprintf('\n');
    fprintf('%4i ',Notavailablelist);
	SV = input('\nin the form of   "[ sv sv . . . ]"  :  ');
% prune away any SV data that is not to be used.
    for jj = Nsats:-1:1
       if ~any(SVIDlist(jj,1) == SV)
          SVIDlist(jj,:) = [];
          pseudoDatahist(:,jj) = [];
          DopplerDatahist(:,jj) = [];
          iPRNflagshist(:,jj) = [];
       end
    end
    clear jj SV
% create the truncated ephem array.
    Nsats = size(SVIDlist,1);
    ephem = zeros(Nsats,24);
    for jj = 1:Nsats
       idumjj = find(SVIDlist(jj,1) == ephemData(:,1));
       ephem(jj,:) = ephemData(idumjj,:);
    end
    clear jj idumjj ephemData
% quit with an error if there are fewer than 4 satellites.
    if Nsats < 4
       disp(' ')
       disp('Error in navvelsolnhist.m: Fewer than 4 satellites')
       disp(' available for navigation and velocity determination.')
       disp(' ')
       return
    end
% set up output arrays
    ECEFposhist = zeros(N,3);
    ECEFvelhist = zeros(N,3);
    latlongalthist = zeros(N,3);
    recCOhist = zeros(N,1);
    recCOdothist = zeros(N,1);
    ttruehist = zeros(N,1);
    GDOPhist = zeros(N,1);
    PDOPhist = zeros(N,1);
    TDOPhist = zeros(N,1);
    HDOPhist = zeros(N,1);
    VDOPhist = zeros(N,1);
    elevhist = zeros(N,Nsats);
    azimhist = zeros(N,Nsats);
    SVIDusedflaghist = zeros(N,Nsats);
    sigmaPRhist = zeros(N,1);
    sigmaDopphist = zeros(N,1);
    Nsatshist = zeros(N,1);
% work through the N samples, solve for the position and
% velocity at each sample, and store the results.
    for k = 1:N
       pseudoRk = pseudoDatahist(k,:)';
       Doppshiftk = DopplerDatahist(k,:)';
       gpsTimek = gpsTimehist(k,1);
       iPRNflagsk = iPRNflagshist(k,:)';
       idumveck = find(iPRNflagsk == 1);
       [posvelOBSk,DOPk,el_azk,SVsusedk,sigmaPRk,sigmaDoppk] = ...
                   solveposvelod(ephem(idumveck,:),...
                                 pseudoRk(idumveck,1),...
                                 Doppshiftk(idumveck,1),guess,...
                                 gpsTimek,ionParam,iflagion,elevmask,...
                                 p,TdegK,hrel,iflagna);
       ECEFposhist(k,:) = posvelOBSk(1,2:4);
       ECEFvelhist(k,:) = posvelOBSk(1,6:8);
       latlongalthist(k,:) = latlong(posvelOBSk(1,2:4));
       recCOhist(k,1) = posvelOBSk(1,5);
       recCOdothist(k,1) = posvelOBSk(1,9);
       ttruehist(k,1) = gpsTimek - posvelOBSk(1,5);
       GDOPhist(k,1) = DOPk(1,1);
       PDOPhist(k,1) = DOPk(2,1);
       TDOPhist(k,1) = DOPk(3,1);
       HDOPhist(k,1) = DOPk(4,1);
       VDOPhist(k,1) = DOPk(5,1);
       elevhist(k,idumveck) = el_azk(:,3)';
       azimhist(k,idumveck) = el_azk(:,4)';
       for jj = 1:Nsats
          if any(SVIDlist(jj,1) == SVsusedk)
             SVIDusedflaghist(k,jj) = 1;
          end
       end
       sigmaPRhist(k,1) = sigmaPRk;
       sigmaDopphist(k,1) = sigmaDoppk;
       Nsatshist(k,1) = size(SVsusedk,1);
       % use the current solution as the next position guess.
       guess = posvelOBSk(1,2:4)';
    end
% clear miscellaneous intermediate results.
    clear k pseudoRk Doppshiftk gpsTimek jj posvelOBSk DOPk ...
          el_azk SVsusedk sigmaPRk sigmaDoppk idumveck
% display statistics about the number of available SVs
    disp(' ')
    disp(['The minimum number of SVs used is ',int2str(min(Nsatshist)),','])
    disp([' the maximum number is ',int2str(max(Nsatshist)),', and the'])
    disp([' mean number is ',num2str(mean(Nsatshist)),'.'])
% compute the average position and its latitude, longitude and altitude
% or input it, or use a user-supplied "average" position.
    if iflagfiles > 5
  	   fprintf('\nEnter "1" if you want to input a nominal position ');
  	   fprintf('\n about which vertical/east/north variations \n');
	   iflagposavg = input(' are recorded, otherwise enter "0"  :  ');
    else
       iflagposavg = 0;
    end
    if iflagposavg == 1
	   fprintf('\nEnter the nominal location about which vertical/');
	   fprintf('\n east/north variations are to be computed \n');
	   latlongaltavg = ...
             input(' in the form "[ latitude longitude altitude ]"  :  ');
       ECEFposavg = ecef(latlongaltavg);
       ECEFposavg = ECEFposavg(:);
    else
       ECEFposavg = [mean(ECEFposhist(:,1));mean(ECEFposhist(:,2));...
                     mean(ECEFposhist(:,3))];
       latlongaltavg = latlong(ECEFposavg');
    end
% compute the transformation from ECEF coordinates to vertical/east/north
% coordinates at the average location.
    phiavg = latlongaltavg(1,1)*pi/180;
    lambdaavg = latlongaltavg(1,2)*pi/180;
    A_VEN_ECEF = ...
       [cos(phiavg),0,sin(phiavg);0,1,0;-sin(phiavg),0,cos(phiavg)]*...
       [cos(lambdaavg),sin(lambdaavg),0;...
          -sin(lambdaavg),cos(lambdaavg),0;0,0,1];
% compute the pseudorange-based position time history in   
% vertical/east/north coordinates with its center at the average position
    VENposhist = [(ECEFposhist(:,1) - ECEFposavg(1,1)),...
                  (ECEFposhist(:,2) - ECEFposavg(2,1)),...
                  (ECEFposhist(:,3) - ECEFposavg(3,1))]*(A_VEN_ECEF');
% compute the velocity in vertical/east/north coordinates.
    VENvelhist = ECEFvelhist*(A_VEN_ECEF');
% compute the velocity-smoothed position time history in
% vertical/east/north coordinates with its center at the average 
% position
    VENposvelsmthhist = integrate_vel_ven(ttruehist,VENvelhist);
    for jj = 1:3
       VENposvelsmthhist_jj = VENposvelsmthhist(:,jj);
        VENposvelsmthhist(:,jj) = VENposvelsmthhist_jj - ...
                            mean(VENposvelsmthhist_jj - VENposhist(:,jj));
    end
    clear jj VENposvelsmthhist_jj
% compute the receiver-clock-error-rate-smoothed receiver clock
% error time history.
    dumvec = 0.5*diff(ttruehist).*...
             (recCOdothist(1:(N-1),1) + recCOdothist(2:N,1));
    recCOsmthhist = cumsum([0;dumvec]);
    clear dumvec
    recCOsmthhist = recCOsmthhist + mean(recCOhist - recCOsmthhist);
% plot the pseudorange-based and velocity-smoothed North vs. East
% positions.
    figure(1)
    hold off
    plot(VENposhist(:,2),VENposhist(:,3),'b-','LineWidth',2)
    hold on
    plot(VENposvelsmthhist(:,2),VENposvelsmthhist(:,3),'r-.','LineWidth',2)
    plot(0,0,'g*','LineWidth',2,'MarkerSize',10)
    hold off
    grid
    axis('equal')
    xlabel('East Position (meters)')
    ylabel('North Position (meters)')
    title(['Northing vs. Easting from * at Lat: ',...
           num2str(latlongaltavg(1,1),'%11.6f'),...
           ' deg, Long: ',num2str(latlongaltavg(1,2),'%11.6f'),...
           ' deg, Alt: ',num2str(latlongaltavg(1,3),'%9.3f'),' m'])
    legend('Pseudorange Solution','Integrated Velocity Solution')
% plot the pseudorange-based and velocity-smoothed vertical displacement
% time histories and receiver clock error time histories.
    figure(2)
    subplot(211)
    hold off
    plot(ttruehist-ttruehist(1,1),VENposhist(:,1),'b-','LineWidth',2)
    hold on
    plot(ttruehist-ttruehist(1,1),VENposvelsmthhist(:,1),'r-.','LineWidth',2)
    hold off
    grid
    ylabel('Vertical Displacement (meters)')
    legend('Pseudorange Solution','Integrated Velocity Solution')
    title(['Vertical Motion Profile Relative to Alt of Nominal Point: ',...
            num2str(latlongaltavg(1,3),'%9.3f'),' m'])
    subplot(212)
    hold off
    plot(ttruehist-ttruehist(1,1),recCOhist,'b-','LineWidth',2)
    hold on
    plot(ttruehist-ttruehist(1,1),recCOsmthhist,'r-.','LineWidth',2)
    hold off
    grid
    xlabel('Elapsed Time (sec)')
    ylabel('Receiver Clock Error (seconds)')
    legend('Pseudorange Solution','Integrated Velocity Solution')
    title('Receiver Clock Error Time History')
% plot the vertical/east/north velocity time history from the velocity
% solution and compute and plot the finite-differenced position solution.
    figure(3)
    tVENvelFDhist = 0.5*(ttruehist(1:(N-1),1) + ttruehist(2:N,1));
    VENvelFDhist = diff(VENposhist)./(diff(ttruehist)*ones(1,3));
    subplot(311)
    hold off
    plot(tVENvelFDhist-ttruehist(1,1),VENvelFDhist(:,1),'b-','LineWidth',2)
    hold on
    plot(ttruehist-ttruehist(1,1),VENvelhist(:,1),'r-.','LineWidth',2)
    hold off
    grid
    ylabel('Vertical Velocity (m/sec)')
    title('Vertical Velocity Time History')
    legend('Finite-Differenced Pseudorange Solution',...
           'Velocity Solution from Doppler')
    subplot(312)
    hold off
    plot(tVENvelFDhist-ttruehist(1,1),VENvelFDhist(:,2),'b-','LineWidth',2)
    hold on
    plot(ttruehist-ttruehist(1,1),VENvelhist(:,2),'r-.','LineWidth',2)
    hold off
    grid
    ylabel('East Velocity (m/sec)')
    title('East Velocity Time History')
    legend('Finite-Differenced Pseudorange Solution',...
           'Velocity Solution from Doppler')
    subplot(313)
    hold off
    plot(tVENvelFDhist-ttruehist(1,1),VENvelFDhist(:,3),'b-','LineWidth',2)
    hold on
    plot(ttruehist-ttruehist(1,1),VENvelhist(:,3),'r-.','LineWidth',2)
    hold off
    grid
    xlabel('Elapsed Time (sec)')
    ylabel('North Velocity (m/sec)')
    title('North Velocity Time History')
    legend('Finite-Differenced Pseudorange Solution',...
           'Velocity Solution from Doppler')       
% plot the DOP time histories. the pseudorange RMS residuals
% time history, and the carrier Doppler shift RMS residuals time history.
    figure(4)
    subplot(311)
    hold off
    plot(ttruehist-ttruehist(1,1),GDOPhist,'b-','LineWidth',2)
    hold on
    plot(ttruehist-ttruehist(1,1),HDOPhist,'g--','LineWidth',2)
    plot(ttruehist-ttruehist(1,1),VDOPhist,'r-.','LineWidth',2)
    hold off
    grid
    ylabel('Dilution of Precision')
    legend('GDOP','HDOP','VDOP')
    title('DOP Time Histories for Position Solution')
    subplot(312)
    hold off
    idum = find(~isnan(sigmaPRhist));
    plot(ttruehist(idum,1)-ttruehist(1,1),sigmaPRhist(idum,1),...
         'b-','LineWidth',2)
    grid
    ylabel('sigmaPR (meters)')
    title('RMS Residual Pseudorange Error Time History') 
    subplot(313)
    hold off
    idum = find(~isnan(sigmaDopphist));
    plot(ttruehist(idum,1)-ttruehist(1,1),sigmaDopphist(idum,1),...
         'b-','LineWidth',2)   
    grid
    xlabel('Elapsed Time (sec)')
    ylabel('sigmaDopp (Hz)')
    title('RMS Residual Doppler Error Time History')    