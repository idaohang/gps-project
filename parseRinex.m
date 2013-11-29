%**************************************************************************
% parseRinex.m
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% This function parses RINEX files into space-delimited arrays appropriate
% for use with GPS course code.
%
%*************************************************************************
function parseRinex()
type = 0;
while(type ~= 1 && type ~= 2)
  fprintf('\nPossible data type to extract the observation file: ');
  fprintf('\nType "1": OBS data');
  fprintf('\nType "2": NAV data');
  type = input('\nEnter data type in the form of "1" or "2": ');
  if(type ~= 1 && type ~= 2)
    fprintf('\nData type is not valid\n\n');
  end
end

if(type == 1)
  rinexObs = input('Observation file name: ', 's');
  fid2 = fopen(rinexObs);
  if (fid2 < 0)
    fprintf('\nRINEX2 observation file not found.\n\n');
    return;
  end
  
  [obs, obs_type] = getobs(fid2);
  % this 'if' block is here for compatibility with 4150 code and the 
  % expectation that L1 C/A pseudoranges are stored in the file 'obs.asc'.
  % If that is not your goal, use only the 'else' block.
  if strcmp(obs_type, 'C1')
    filename = 'obs.asc';
    save 'obs.asc' obs -ascii -double
  else
    idx = find(rinexObs == '.');
    if(isempty(idx))
      file = rinexObs;
    else
      file = rinexObs(1:idx-1);
    end
    filename = sprintf('%s_%s.mat',file,obs_type);
    save(filename,'obs');
  end
  fclose(fid2);
  
  fprintf(sprintf('\n%s has been created successfully.\n\n',filename));
elseif(type == 2)
  fprintf('Enter RINEX2 navigation and observation files.\n');
  rinexNav = input('Navigation file name: ', 's');
  fid = fopen(rinexNav);
  if (fid < 0)
    fprintf('\nRINEX2 navigation file not found.\n\n');
    return;
  end
  [E, I] = getephem(fid, 0);
  save 'ion.asc' I -ascii -double;
  save 'ephem.asc' E -ascii -double;
  fclose(fid);
  fprintf('\nephem.asc and ion.asc created successfully\n\n');
end;
return


%**************************************************************************
% getephem.m
%
% this code segment opens RINEX2 navigation files and
% generates an ephem matrix that is compatible with
% the gps course code (ECE 415)
%
% input:    'infile' the file id (fid) of a valid RINEX2 navigation file
%            use fopen to generate a fid
%           'sats' the list of satellites recoqnized in the observation
%           file
% output:   'E' an ephem matrix
%
%**************************************************************************

function [E, I] = getephem(infile, sats)
if(sats == 0)
  sats = [1:32];
end
E = zeros(32,24);

% get the header info
buff=getline(infile);
IONA = [0; 0; 0; 0];
IONB = [0; 0; 0; 0];
while (~isempty(buff))
  if ~isempty(strfind(buff,'ION ALPHA'))
    buff = strrep(buff, 'D', 'E');
    IONA = sscanf(buff,'   %e%e%e%e',82);
    
    buff = getline(infile);
    buff = strrep(buff, 'D', 'E');
    IONB = sscanf(buff,'   %e%e%e%e',82);
  end;
  if (~isempty(strfind(buff,'END OF HEADER')))
    break;
  end
  buff=getline(infile);
end

I = [IONA IONB];

% check for end of file
if (feof(infile))
  return;
end

while (~feof(infile))
  buff = getline(infile);
  % add the satellite to the appropriate row
  % this row has info about the time and date
  if (feof(infile))
    break;
  end
  buff = strrep(buff, 'D', 'E');
  ROW1 = sscanf(buff,'%d %d %d %d %d %d %e%e%e%e');
  
  % get the gpstime
  gpstime = timeConversion(ROW1(2),ROW1(3),ROW1(4),ROW1(5),ROW1(6),ROW1(7));
  
  % Treat ephemeris as valid 2 hrs before or 4 hrs after
  % time of ephemeris.
  buff = getline(infile); buff = strrep(buff, 'D', 'E'); ROW2 = sscanf(buff,'   %e%e%e%e',82);
  buff = getline(infile); buff = strrep(buff, 'D', 'E'); ROW3 = sscanf(buff,'   %e%e%e%e',82);
  buff = getline(infile); buff = strrep(buff, 'D', 'E'); ROW4 = sscanf(buff,'   %e%e%e%e',82);
  buff = getline(infile); buff = strrep(buff, 'D', 'E'); ROW5 = sscanf(buff,'   %e%e%e%e',82);
  buff = getline(infile); buff = strrep(buff, 'D', 'E'); ROW6 = sscanf(buff,'   %e%e%e%e',82);
  buff = getline(infile); buff = strrep(buff, 'D', 'E'); ROW7 = sscanf(buff,'   %e%e%e%e',82);
  buff = getline(infile); buff = strrep(buff, 'D', 'E'); ROW8 = sscanf(buff,'   %e%e%e%e',82);
  
  % fill in a row of E
  sv = ROW1(1);
  
  % dont extract ephemerides for satelites not listed in the
  % observation file.
  if(~isempty(find(sats==sv, 1)))
    E(sv,1) = sv;
    % get the gpsweek
    [~, E(sv,2)] = timeConversion(ROW1(2),ROW1(3),ROW1(4),0,0,0);  %need to convert a GMT time to GPS
    E(sv,4) = ROW4(1); %toe
    E(sv,3) = E(sv,2)*604800 + E(sv,4);  %weeks*604800 + toe
    E(sv,5) = ROW3(2); %ecc;
    E(sv,6) = ROW3(4); %sqrtA;
    E(sv,7) = ROW5(3); %w0
    E(sv,8) = ROW2(4);  %M0
    E(sv,9) = ROW4(3); %OMEGA
    E(sv,10) = ROW5(4); %OMEGA dot
    E(sv,11) = ROW2(3);  %deltaN
    E(sv,12) = ROW5(1); %incl;
    E(sv,13) = ROW6(1); %idot;
    E(sv,14) = ROW3(1); %Cuc;
    E(sv,15) = ROW3(3); %Cus;
    E(sv,16) = ROW5(2); %Crc
    E(sv,17) = ROW2(2); %Crs;
    E(sv,18) = ROW4(2); %Cic;
    E(sv,19) = ROW4(4); %Cis;
    E(sv,20) = ROW1(8); %af0
    E(sv,21) = ROW1(9); %af1
    E(sv,22) = ROW1(10); %af2
    E(sv,23) = ROW7(3); %tgd
    E(sv,24) = gpstime; %toc
  end
end
% now, remove the zeros
idx = E(:,1) ~= 0;
E = E(idx,:);

frewind(infile);

return;


%***************************************************************************
%* Get the next line from the input file (up to 80 characters).
% This SKIPS ALL COMMENT LINES
%**************************************************************************
function line = getline(file)
line = fgets(file);
while ~isempty(strfind(line,'COMMENT'))
  line = fgets(file);
end
return;

%**************************************************************************
function [gpstime,gpsweek] = timeConversion(year,month,day,hour,minute,second)

%convert day from year, month, day to Julian day */
%epoch for Julian day is Jan 1, 4713 BC */
%handle different leap years between Julian and Gregorian calenders */
if ((month == 1) || (month == 2))
  year = year - 1;
  month = month + 12;
end
%calculate Julian date */
if(year < 70)
  y = 19.0 + 1;
else
  y = 19.0;
end

julianDate = (2.0-y+floor(y/4.0)) + floor(365.25*(y*100.0+year))  + ...
  floor(30.6001*(month+1.0)) + day + 1720994.5;

% determine day of week
a = floor(julianDate + 0.5);
dayOfWeek = floor(a) - 7.0 * floor(a/7.0) + 1.0;
% determine gps time
gpstime = dayOfWeek * 86400.0 + hour * 3600.0 + minute * 60.0 + second;

gpsweek = floor((julianDate - 2444244.5)/7.0);

return;

%**************************************************************************
% getobs.m
%
% this code segment opens RINEX2 observation files and
% generates an ephem matrix that is compatible with
% the gps course code (ECE 415)
% input:    'infile' the file id (fid) of a valid RINEX2 navigation file
%            use fopen to generate a fid
%
%
% output:   'obs' matrix
%**************************************************************************

function [obs, obs_type] = getobs(infile)
format compact;

% get the header info
buff=getline(infile);

numobs = 0;

while (~isempty(buff))
  if (~isempty(strfind(buff,'END OF HEADER')))
    break;
  elseif(~isempty(strfind(buff,'INTERVAL')))
    interval = round(sscanf(buff,'%f'));
  elseif(~isempty(strfind(buff,'TYPES OF OBSERV')))
    numobs = sscanf(buff,'%d');
    obstype_array  = ['L1'; 'L2'; 'P1'; 'P2'; 'C1'; 'C2'; 'S1'; 'S2'; 'D1'; 'D2'];
    obstype = zeros(10,1);
    obsorder = zeros(10,1);
    for x=1:length(obstype_array)
      temp = strfind(buff,obstype_array(x,:));
      if(~isempty(temp))
        obstype(x,1) = 1;
        obsorder(x,1) = temp;
      end
    end
  end
  buff=getline(infile);
end

MAXOBS = 30;  %maximum of 14 channels of data
if(~exist('interval','var'))
  interval = 30; %assume modulo-30 interval
end
obs = zeros(ceil(86400/interval),MAXOBS);

%ask the user for what data they want
idx = find(obstype==1);
idx0 = find(obstype==0);
if(isempty(idx))
  fprintf('\nI could not interpret the RINEX file\n\n');
  return;
end
buff = [];
for x=1:length(idx)
  buff = [buff sprintf('%d: %s   ',x,obstype_array(idx(x),:))];
end
fprintf(sprintf('\nThe following data types are available from the obs file:\n%s\n',buff));
fprintf(' L1, L2 are phase\n C1, C2 are pseudorange\n D1, D2 are doppler\n S1, S2 are signal strengths.\n');
temp = inf;
while(temp>length(idx));
  temp = input('Please select one of the file types: ');
end
obs_type = obstype_array(idx(temp),1:2);
[y,~] = sort(obsorder);
obstype = find(y==obsorder(idx(temp)))-length(idx0);

% check for end of file
if (feof(infile))
  return;
end

k = 0;
sats_arr = zeros(1,32);
while (~feof(infile))
  k = k + 1;
  flag = 0;
  buff = getline(infile);

  % add the satellite to the appropriate row
  % this row has info about the time and date
  if(feof(infile))
    break;
  end
  idx = find(buff=='G', 1);
  if(isempty(idx))
    ROW1 = sscanf(buff,'%d %d %d %d %d %e %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d');
    flag = 1;
  else
    ROW1 = sscanf(buff,'%d %d %d %d %d %e %d %dG%dG%dG%dG%dG%dG%dG%dG%dG%dG%dG%dG%dG%dG%d');
  end
  % check for event flag missing epoch
  if numel(ROW1)<3 %found an event flag in the file 
    for jjj = 1:ROW1(2)
      buff = fgets(infile);
      fprintf('FOUND EVENT: %s', buff);
    end
    continue
  end
  numsats = ROW1(8);
  if(numsats>12)
    sats_temp = ROW1(9:9+12-1)';
    buff = getline(infile);
    if(flag)
      ROW2 = sscanf(buff,'                                 %d %d %d');
    else
      ROW2 = sscanf(buff,'                                G%dG%dG%d');
    end
    sats_temp = [sats_temp ROW2(1:numsats-12)'];
  else
    sats_temp = ROW1(9:9+numsats-1)';
  end
  
  %added more intelligent sats extraction algorithm
  
  for iii=1:length(sats_temp)
    if(sats_temp(iii)~=0) %for some receivers, PRN 0 indicates bad data or dummy SV
      sats_arr(sats_temp(iii)) = 1;
    end
  end;
  
  % get the gpstime
  gpstime = mod(timeConversion(ROW1(1),ROW1(2),ROW1(3),ROW1(4),ROW1(5),ROW1(6)),604800);
  
  
  sample = zeros(1,30);
  sample(1:2) = [k gpstime];
  obssamples = nan(numsats,1);
  for i=1:numsats
    observations = nan(numobs,1);
    count = 0;
    buff = getline(infile);
    
    % rinex data observations are in 14.3f format, with the possiblity of
    % 1 digit loss-of-lock & signal strength indicators following each obs,
    % so actual observables are 16 chars apart and 14 chars wide.
    for x=1:16:length(buff)-14
      count = count+1;
      temp = sscanf(buff(x:x+13), '%f', 1);
      if(isempty(temp))
        temp = nan;
      end
      observations(count) = temp;
    end
    if(numobs > 5)
      buff = getline(infile);
      for x=1:16:length(buff)-14
        count = count+1;
        temp = sscanf(buff(x:x+13), '%f', 1);
        if(isempty(temp))
          temp = nan;
        end
        observations(count) = temp;
      end
    end
    obssamples(i) = observations(obstype);
  end
  
  % build up the sample array
  sample(3:2:2+2*length(sats_temp)) = sats_temp;
  sample(4:2:2+2*length(sats_temp)+1) = obssamples;
  
  obs(k,:) = sample;
  
end
obs = obs(1:k,:);
frewind(infile);

return;



