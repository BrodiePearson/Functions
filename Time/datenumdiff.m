function [d] = datenumdiff(d1,d2,unit)
% calculates time difference between two datenumbers, d1 and d2
% d1 = initial time, d2 = final time
% type specifies the unit returned
% unit = 'hour', 'second', or 'minute' or 'string'

d = datetime(d2, 'ConvertFrom', 'datenum') - datetime(d1, 'ConvertFrom', 'datenum');

if strcmp(unit, 'day') == 1
    d = days(d);
elseif strcmp(unit, 'hour') == 1
    d = hours(d);
elseif strcmp(unit, 'minute') == 1
    d = minutes(d);
elseif strcmp(unit, 'second') == 1
    d = seconds(d);   
end