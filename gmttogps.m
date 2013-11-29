% gmttogps.m    (actual file name: gmttogps.m)
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% this GPS utility converts a GMT time and date into a GPS time
% of week (seconds into the week from midnight on Sunday)
%
% input: 'date' which contains the GMT date in a string format
%       mmddyy
%           'gmtTime' which contains the GMT time in a string format
%       hhmmss
%
% output: 'gptTime' which contains the converted GPS time of week
%                               [ gpsTime ]
%
function gpsTime = gmttogps(date,gmtTime)
% define physical constants
    constant;
% convert the strings 'date' and 'gmtTime' into the components 
% year, month, day, hour, minutes, and seconds
    year = str2double(date(5:6)) + 2000;
    day = str2double(date(3:4));
    month = str2double(date(1:2));
    second = (gmtTime(5:end));
    temp = 0;
    if(length(second)>=4)
        second = [second(1:2), second(4:end)];
        power = fliplr(10.^((0:1:length(second))-(length(second)-1)));
        for x=1:length(second)
            temp = temp+(second(x)-48)*power(x);
        end;
        second = temp;
    else
        second = str2double(gmtTime(5:6));
    end
    minute = str2double(gmtTime(3:4));
    hour = str2double(gmtTime(1:2));
% convert day from year, month, day to Julian day number
    % the epoch for the Julian day is taken as Greenwich mean noon
    % of January 1st, 4713 BC
    % handle different leap years between Julian and Gregorian 
    % calendars
    if (month == 1 || month == 2)
        year = year - 1;
        month = month + 12;
    end
    % assume date to be converted occurs after Julian versus
    % Gregorian calendar which occurred 15 Oct 1582; calculate 
    % Julian date 'jdn'
    y = fix(year / 100);
    jdn = (2 - y + fix(y / 4)) + fix(365.25 * year) + ...
        fix(30.6001 * (month + 1)) + day + 1720994.5;
% determine day of week starting with N = 0 for Sunday
    a = fix(jdn + 0.5);
    N = mod(a - 7 * fix(a / 7) + 1,7);
% compute the number of seconds from Sunday midnight GMT; return 
% the calculated GPS time 
    gpsTime = N*24*60*60 + hour*60*60 + minute*60 + second;
    % compensate for leap second differential between GMT time
    % and GPS time
    gpsTime = gpsTime + 16;
% return 'gpsTime'
    return;
