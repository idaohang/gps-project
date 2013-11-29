% pseudorangeerrors2.m		(actual file name: pseudorangeerrors2.m)
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% < MULTI-SATELLITE OVERDETERMINED SOLUTION FOR CLOCK ERROR ONLY
%   FOLLOWED BY CALCULATION OF PSEUDORANGE ERRORS INCLUDING
%   IONOSPHERE AND NEUTRAL ATMOSPHERE EFFECTS.>
%
% This program determines a receiver clock error solution directly
% from raw pseudo-ranges or from carrier-smoothed pseudoranges,
% satellite ephemerides, and a known receiver position. The user is 
% required to supply the known receiver location, a list of
% satellite SV ids to use in the clock error solution, and the
% pressure, temperature, and relative humidity at the receiver.
%
% The files required for pseudorangeerrors2.m are:
%
%		* constant.m
%       * deltl1ionocalc.m
%       * deltnacalc.m
%		* ecef.m
%       * elevazim.m
%		* findsatclock.m
%       * latlong.m
%		* solveclockod2.m
%		
% Ephemeris data and pseudo-range data are needed for calculating
% the corrected pseudo-ranges, the satellite positions, and the
% final clock error solution. The program parseRinex.m must have 
% already been used to format data into two files; ephem.asc and 
% obs.asc. If L1 beat carrier phase data are to be used too, then
% a third call of parseRinex.m must have formated data in to a file
% with a name of the form: 'XXXXX_L1.mat', where XXXXX is the
% name of the original RINEX2 observation data file from
% which parseRinex.m retrieved the beat carrier phase data.
%
% This script also computes the azimuths and elevations of the 
% satellites.
%
% pseudorangeerrors2.m consists of the following steps and function calls:
%
%		* opens data files - ephem.asc, ion.asc, and obs.asc
%       * inputs neutral atmosphere pressure (p), temparature (TdegC),
%         and relative humidity (hrel) data and computes temperature
%         in degrees Kelvin (TdegK)
%       * inputs flags that determine whether ionosphere delay
%         corrections and neutral atmosphere delay corrections
%         are applied, during both solution for clock errors
%         and calculation of pseudorange errors.
%		* inputs decision whether L1 beat carrier phase data to be used
%       * inputs name of beat carrier phase data file and opens it
%         and translates it into equivalent accumulated delta ranges.
%		* inputs the observation station location
%       * sorts through the observables to determine which satellites
%           have pseudorange data available for the entire span
%           formats pseudorange data and beat carrier phase data
%           into arrays for later processing and sets up 
%           indexing arrays controlling access to the various 
%           available satellites ephemeris, pseudo-range, 
%           and beat carrier-phase data.
%		* inputs satellite SVs to be used for clock error
%           solution
%		* Removes mean difference between accumulated delta range and
%           pseudo-range time histories to turn accumulated delta 
%           range into carrier-smoothed code phase if beat carrier 
%           phase is considered.
%		* cycles through all samples, calling solveclockod2.m,
%           which iteratively calculates the receiver's clock
%			error from satellite locations, pseudo-ranges
%           or carrier-smoothed pseudo-ranges, and observer
%           location, computes pseudo-range residual errors,
%           azimuths, elevations, ionospheric delays, and
%           neutral atmosphere delays, and stores residual 
%           error results, satellite clock error results,
%           azimuths, elevations, ionospheric delays, and
%           neutral atmosphere delay in appropriate arrays.
%       * re-iterates the last 2 bullet points if beat carrier phase
%           data have been used and if ionospheric corrections are 
%           used.  This re-iteration is done in order to remove 
%           twice the mean ionospheric correction during the 
%           carrier phase bias calculation.
%		* computes RMS pseudorange errors and saves results in named 
%         arrays as follows:
%
%           'gpsTimehist'      N-by-1 vector of receiver clock times
%                              (sec) of samples.
%           'recCOhist'        N-by-1 vector of receiver clock errors
%                              (sec) of samples.  Note that the
%                              true GPS times of the samples are
%                              gpsTimehist - recCOhist.
%           'SVIDusedlist'     The list of PRN numbers that have been
%                              used to compute the receiver clock error.
%           'obsPos'           3-by-1 ECEF Cartesian position vector
%                              of receiver (meters) as determined
%                              from user input obsPoslatlongalt
%           'obsPoslatlongalt' 1-by-3 vector of WGS-84 latitude (deg)
%                              longitude (deg) and altitude (meters)
%                              as input by user.
%           'NPRNXX'           The number of samples for which data have
%                              been processed for PRNXX.
%           'isampsPRNXX'      The NPRNXX-by-1 vector of time sample 
%                              indices at whith the pseudo-range  
%                              errors for PRNXX in prerrorPRNXX and
%                              the carrier phase errors for PRNXX
%                              in cperrorPRNXX apply.
%           'prerrorPRNXX'     NPRNXX-by-1 vector of pseudo-range
%                              errors (meters) for PRNXX at the sample
%                              times gpsTimehist(isampsPRNXX,1).
%           'cperrorPRNXX'     NPRNXX-by-1 vector of beat carrier phase
%                              errors (meters) for PRNXX at the sample
%                              times gpsTimehist(isampsPRNXX,1).  These
%                              are actually errors in the carrier-
%                              smoothed pseudorange.
%           'azimuthPRNXX'     NPRNXX-by-1 vector of azimuths
%                              (deg) for PRNXX at the sample
%                              times gpsTimehist(isampsPRNXX,1).
%           'elevationPRNXX'   NPRNXX-by-1 vector of elevations
%                              (deg) for PRNXX at the sample
%                              times gpsTimehist(isampsPRNXX,1).
%           'deltl1ionoPRNXX'  NPRNXX-by-1 vector of modeled ionospheric 
%                              delays (sec) for PRNXX at the sample
%                              times gpsTimehist(isampsPRNXX,1).
%           'deltnaPRNXX'      NPRNXX-by-1 vector of modeled neutral 
%                              atmosphere delays (sec) for PRNXX 
%                              at the sample times 
%                              gpsTimehist(isampsPRNXX,1).
%
%           Note: outputs with the names NPRNXX, isampsPRNXX, 
%           prerrorPRNXX, cperrorPRNXX, azimuthPRNXX, elevationPRNXX, 
%           deltl1ionoPRNXX, and deltnaPRNXX will be repeated for 
%           all available PRN numbers in the data.

%
% < MULTI-SATELLITE OVERDETERMINED SOLUTION FOR CLOCK ERROR ONLY
%   FOLLOWED BY CALCULATION OF PSEUDORANGE ERRORS INCLUDING
%   IONOSPHERE AND NEUTRAL ATMOSPHERE EFFECTS.>
%

% clear Matlab workspace
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
% determine number of samples and the receiver clock GPS times
    N = size(pseudoData,1);
    gpsTimehist = pseudoData(:,2);
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
% input the decision of whether to include beat carrier phase data
	fprintf('\nEnter a "1" if beat carrier phase data are to \n');    
	iflagcpdata = input('be used.  Enter a "0" otherwise : ');
% retrieve the beat carrier phase data and translate it into
% meters.  put the raw pseudo-range data in pseudoData2 and
% the range-equivalent beat carrier phase data in pseudoData.
    if iflagcpdata == 1
       cpdatafilename = ...
            input('\nEnter name of beat carrier phase .mat file : ','s');
       load(cpdatafilename);
       pseudoData2 = pseudoData;
       pseudoData = obs;
       clear obs
       if size(pseudoData,1) ~= N
          error(['Error in pseudorangeerrors2.m:  The pseudo-range',...
                 ' observables file and the beat carrier phase',...
                 ' observables do not have the same number of samples.'])
       end
       constant
       pseudoData(:,4:2:end) = pseudoData(:,4:2:end)*lambdaL1;
       clear AA OmegaE degrad f0 lambdaL1 muearth BB c esquare ...
             fL1 leapSeconds
       if norm(pseudoData(:,1:2) - pseudoData2(:,1:2)) > 0
          error(['Error in pseudorangeerrors2.m:  The pseudo-range',...
                 ' observables file and the beat carrier phase',...
                 ' observables file have different sample numbers',...
                 ' or different sample times.'])
       end
       if norm(pseudoData(:,3:2:end) - pseudoData2(:,3:2:end)) > 0
          error(['Error in pseudorangeerrors2.m:  The pseudo-range',...
                 ' observables file and the beat carrier phase',...
                 ' observables file have different PRN numbers.'])
       end
       clear cpdatafilename
    else
       pseudoData2 = [];
    end
% input receiver's position in latitude-longitude
% coordinates 
	fprintf('\nEnter the precise location of the observation station \n');
	obsPoslatlongalt = ...
        input(['	in the form "[ latitude(deg) longitude(deg)',...
               ' altitude(m) ]"  :  ']);
    obsPos = ecef(obsPoslatlongalt);
    obsPos = obsPos(:);
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
    clear nsatsdum jj
% create pure pseudorange data arrays and indexing arrays for the
% available pseudorange data.  Pseudo-range and carrier-smoothed
% pseudorange data in the arrays pseudoDatahist and
% pseudoData2hist are valid only if the corresponding
% entries of iPRNflagshist contain 1 values.  Each
% column corresponds to the GPS satellite whose SVID is
% contained in the corresponding row of SVIDlist, and each
% row corresponds to the time in the corresponding row
% of gpsTimehist.
    Nsats = size(SVIDlist,1);
    iPRNflagshist = zeros(N,Nsats);
    pseudoDatahist = zeros(N,Nsats);
    if iflagcpdata == 1
       pseudoData2hist = zeros(N,Nsats);
    end
    for k = 1:N
       SVIDlistk = pseudoData(k,3:2:end)';
       for jj = 1:Nsats
          idumkjj = find(SVIDlistk == SVIDlist(jj,1));
          if size(idumkjj,1) == 1
             iPRNflagshist(k,jj) = 1;
             idumkjj_mod = 2*(idumkjj + 1);
             pseudoDatahist(k,jj) = pseudoData(k,idumkjj_mod);
             if iflagcpdata == 1
                pseudoData2hist(k,jj) = pseudoData2(k,idumkjj_mod);
             end
          end
       end
    end
    clear k jj idumkjj idumkjj_mod pseudoData pseudoData2
% create a short list of satellites whose pseudorange data are
% available at all sample times.
    SVIDshortlist = SVIDlist;
    for jj = Nsats:-1:1
       if any(iPRNflagshist(:,jj) == 0)
          SVIDshortlist(jj,:) = [];
       end
    end
% determine which satellites will be used in the receiver clock error
% solution
	fprintf('\nEnter the satellites to be used in the clock error solution ');
	fprintf('\nEnter at least one:');
	fprintf('\nfrom the list of SV ids --> ');
    fprintf('%d ',SVIDshortlist);
    fprintf('\n');
	SVIDusedlist = input('	in the form of   "[ sv sv . . . ]"  :  ');
    SVIDusedlist = SVIDusedlist(:);
    if size(SVIDusedlist,1) < 1
       error(['Error in pseudorangeerrors2.m:  Must use at least',...
              ' one SV in order to determine receiver clock error.'])
    end
% remove the mean error between the raw pseudorange and the accumulated
% delta range in order to synthesize, in effect, carrier-smoothed
% pseudorange. 
    if iflagcpdata == 1
       for jj = 1:Nsats
          isampsjj = find(iPRNflagshist(:,jj) == 1);
          pseudoDatahistjj = pseudoDatahist(isampsjj,jj);
          pseudoData2histjj = pseudoData2hist(isampsjj,jj);
          pseudoDatahist(isampsjj,jj) = pseudoDatahistjj + ...
                       mean(pseudoData2histjj - pseudoDatahistjj);
       end
       clear jj isampsjj pseudoDatahistjj pseudoData2histjj
    end
% initialize arrays to store results for the N samples.
    recCOhist = zeros(N,1);
    prerrorPRNshist = zeros(N,Nsats);
    itypepR = 1;
    itypepR2 = [];
    if iflagcpdata == 1
       cperrorPRNshist = zeros(N,Nsats);
       itypepR = 2;
       itypepR2 = 1;
    end
    azimuthPRNshist = zeros(N,Nsats);
    elevationPRNshist = zeros(N,Nsats);
    deltl1ionoPRNshist = zeros(N,Nsats);
    deltnaPRNshist = zeros(N,Nsats);
% cycle through the N samples, compute the receiver clock error
% for each sample and the corresponding pseudo-range error
% using solveclockod2.m, and store the results in recCOhist, 
% prerrorPRNshist, cperrorPRNshist.  Store auxiliary
% results in azimuthPRNshist, elevationPRNshist, 
% deltl1ionoPRNshist, and deltnaPRNshist. 
   for k = 1:N
      gpsTimek = gpsTimehist(k,1);
% set up inputs to solveclockod.m using only those SVs that
% have valid data at this time.  Also, set the flags in
% iflagpRk only for those SVs that are to be used for
% computing the receiver clock error.
      iPRNflagsk = iPRNflagshist(k,:)';
      idumSVsk = find(iPRNflagsk == 1);
      SVIDlistk = SVIDlist(idumSVsk,1);
      pseudoRk = pseudoDatahist(k,idumSVsk)';
      if iflagcpdata == 1
         pseudoR2k = pseudoData2hist(k,idumSVsk)';
      else
         pseudoR2k = [];
      end
      Nsatsk = size(SVIDlistk,1);
      ephemk = zeros(Nsatsk,24);
      iflagpRk = zeros(Nsatsk,1);
      for jj = 1:Nsatsk
         SVIDkjj = SVIDlistk(jj,1);
         idumkjj = find(SVIDkjj == ephemData(:,1));
         ephemk(jj,:) = ephemData(idumkjj,:);
         if any(SVIDkjj == SVIDusedlist)
            iflagpRk(jj,1) = 1;
         end
      end
      [recCOk,pseudoerrork,pseudoerror2k,azimveck,elevveck,...
                   deltl1ionoveck,deltnaveck] = ...
           solveclockod2(ephemk,pseudoRk,itypepR,iflagpRk,obsPos,...
                         gpsTimek,ionParam,iflagion,p,TdegK,hrel,...
                         iflagna,pseudoR2k,itypepR2);
% store the results for this sample
      recCOhist(k,1) = recCOk;
      if iflagcpdata == 1
         cperrorPRNshist(k,idumSVsk) = pseudoerrork';
         prerrorPRNshist(k,idumSVsk) = pseudoerror2k';
      else
         prerrorPRNshist(k,idumSVsk) = pseudoerrork';
      end
      azimuthPRNshist(k,idumSVsk) = azimveck';
      elevationPRNshist(k,idumSVsk) = elevveck';
      deltl1ionoPRNshist(k,idumSVsk) = deltl1ionoveck';
      deltnaPRNshist(k,idumSVsk) = deltnaveck';
   end
% re-calculate the carrier-smoothed code phase if it is used and
% if ionospheric delay is used, and re-do the calculations using
% the new re-smoothed code phase.
   if (iflagcpdata == 1) & (iflagion == 1)
% re-do the removal of the mean error between the raw pseudorange and 
% the accumulated delta range in order to synthesize, in effect, 
% carrier-smoothed pseudorange. This re-removal accounts for the
% fact that the modeled ionosphere delay affects pseudorange and
% accumulated delta range with the opposite signs.
       constant
       for jj = 1:Nsats
          isampsjj = find(iPRNflagshist(:,jj) == 1);
          pseudoDatahistjj = pseudoDatahist(isampsjj,jj);
          pseudoData2histjj = pseudoData2hist(isampsjj,jj);
          deltl1ionoPRNshistjj = deltl1ionoPRNshist(isampsjj,jj);
          pseudoDatahist(isampsjj,jj) = pseudoDatahistjj + ...
                       mean(pseudoData2histjj - pseudoDatahistjj ...
                            - (2*c)*deltl1ionoPRNshistjj);
       end
       clear jj isampsjj pseudoDatahistjj pseudoData2histjj ...
             deltl1ionoPRNshistjj
       clear AA OmegaE degrad f0 lambdaL1 muearth BB c esquare ...
             fL1 leapSeconds
% re-do the calculations of the receiver clock error and of the
% pseudorange errors, elevations, azimuths, ionospheric delays,
% and neutral atmosphere delays.
      for k = 1:N
         gpsTimek = gpsTimehist(k,1);
% re-set-up inputs to solveclockod.m using only those SVs that
% have valid data at this time.  Also, set the flags in
% iflagpRk only for those SVs that are to be used for
% computing the receiver clock error.
         iPRNflagsk = iPRNflagshist(k,:)';
         idumSVsk = find(iPRNflagsk == 1);
         SVIDlistk = SVIDlist(idumSVsk,1);
         pseudoRk = pseudoDatahist(k,idumSVsk)';
         if iflagcpdata == 1
            pseudoR2k = pseudoData2hist(k,idumSVsk)';
         else
            pseudoR2k = [];
         end
         Nsatsk = size(SVIDlistk,1);
         ephemk = zeros(Nsatsk,24);
         iflagpRk = zeros(Nsatsk,1);
         for jj = 1:Nsatsk
            SVIDkjj = SVIDlistk(jj,1);
            idumkjj = find(SVIDkjj == ephemData(:,1));
            ephemk(jj,:) = ephemData(idumkjj,:);
            if any(SVIDkjj == SVIDusedlist)
            iflagpRk(jj,1) = 1;
            end
         end
         [recCOk,pseudoerrork,pseudoerror2k,azimveck,elevveck,...
                      deltl1ionoveck,deltnaveck] = ...
              solveclockod2(ephemk,pseudoRk,itypepR,iflagpRk,obsPos,...
                            gpsTimek,ionParam,iflagion,p,TdegK,hrel,...
                            iflagna,pseudoR2k,itypepR2);
% re-store the results for this sample
         recCOhist(k,1) = recCOk;
         if iflagcpdata == 1
            cperrorPRNshist(k,idumSVsk) = pseudoerrork';
            prerrorPRNshist(k,idumSVsk) = pseudoerror2k';
         else
            prerrorPRNshist(k,idumSVsk) = pseudoerrork';
         end
         azimuthPRNshist(k,idumSVsk) = azimveck';
         elevationPRNshist(k,idumSVsk) = elevveck';
         deltl1ionoPRNshist(k,idumSVsk) = deltl1ionoveck';
         deltnaPRNshist(k,idumSVsk) = deltnaveck';
      end
   end
% remove unneeded interim quantities    
   clear k gpsTimek iPRNflagsk idumSVsk SVIDlistk pseudoRk pseudoR2k ...
         Nsatsk ephemk iflagpRk jj SVIDkjj idumkjj recCOk ...
         pseudoerrork pseudoerror2k pseudoDatahist ...
         pseudoData2hist azimveck elevveck deltl1ionoveck deltnaveck
% store the pseudo-range error results for the different
% satellites along with the elevation, azimuth, ionospheric
% delay, and neutral atmosphere delay results.  Also, compute
% pseudorange RMS errors, both for used satellites and for all
% satellites.
   sigmaPRsqtotsum = 0;
   nPRsqtotsum = 0;
   sigmaPRsqusedsum = 0;
   nPRsqusedsum = 0;
   for jj = 1:Nsats
      SVjj = SVIDlist(jj,1);
      isampsjj = find(iPRNflagshist(:,jj) == 1);     
      prerrorjj = prerrorPRNshist(isampsjj,jj);
      Njj = size(isampsjj,1);
      sumprerrorjjsq = sum(prerrorjj.^2);
      sigmaPRsqtotsum = sigmaPRsqtotsum + sumprerrorjjsq;
      nPRsqtotsum = nPRsqtotsum + Njj;
      if any(SVjj == SVIDusedlist)
         sigmaPRsqusedsum = sigmaPRsqusedsum + sumprerrorjjsq;
         nPRsqusedsum = nPRsqusedsum + Njj;
      end
      if SVjj < 10
         SVjjtext = ['0',int2str(SVjj)];
      else
         SVjjtext = int2str(SVjj);
      end
      eval(['NPRN',SVjjtext,' = Njj;']);
      eval(['isampsPRN',SVjjtext,' = isampsjj;']);
      eval(['prerrorPRN',SVjjtext,' = prerrorjj;']);
      if iflagcpdata == 1
         cperrorjj = cperrorPRNshist(isampsjj,jj);
         eval(['cperrorPRN',SVjjtext,' = cperrorjj;']);
      end
      azimuthPRNjj = azimuthPRNshist(isampsjj,jj);
      eval(['azimuthPRN',SVjjtext,' = azimuthPRNjj;']);
      elevationPRNjj = elevationPRNshist(isampsjj,jj);
      eval(['elevationPRN',SVjjtext,' = elevationPRNjj;']);
      deltl1ionoPRNjj = deltl1ionoPRNshist(isampsjj,jj);
      eval(['deltl1ionoPRN',SVjjtext,' = deltl1ionoPRNjj;']);
      deltnaPRNjj = deltnaPRNshist(isampsjj,jj);
      eval(['deltnaPRN',SVjjtext,' = deltnaPRNjj;']);
   end
   clear jj SVjj isampsjj prerrorjj Njj SVjjtext cperrorjj
   clear SVIDlist SVIDshortlist cperrorPRNshist ephemData ...
         iPRNflagshist iflagcpdata prerrorPRNshist azimuthPRNshist ...
         elevationPRNshist azimuthPRNjj elevationPRNjj ...
         deltl1ionoPRNshist deltl1ionoPRNjj deltnaPRNshist ...
         deltnaPRNjj sumprerrorjjsq
% compute the RMS pseudorange errors and display them.  the subtraction
% of N in the denominator comes from the fact that N independent
% receiver clock errors have been calculated in the process of
% determining the pseudorange errors.
   sigmaPRused = sqrt(sigmaPRsqusedsum/(nPRsqusedsum - N));
   sigmaPRtot = sqrt(sigmaPRsqtotsum/(nPRsqtotsum - N));
   clear sigmaPRsqtotsum nPRsqtotsum sigmaPRsqusedsum nPRsqusedsum
   fprintf(['\nRMS pseudorange error for satellites used in clock',...
            ' solution: sigmaPRused = ',num2str(sigmaPRused),' meters'])
   fprintf(['\nRMS pseudorange error for all satellites:',...
            ' sigmaPRtot = ',num2str(sigmaPRtot),' meters \n'])
% tell operator that end has occurred.
   fprintf(['\nData for ',int2str(Nsats),' GPS satellites have been'])
   fprintf(['\n processed over ',int2str(N),' samples.\n'])
   clear N Nsats