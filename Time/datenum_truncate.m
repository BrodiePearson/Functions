function datenum_truncated = datenum_truncate(datenum_in, time_unit_string)

[yrs, mons, days, hrs, mins, secs] = datevec(datenum_in);

if strcmp(time_unit_string,'year')
    datenum_truncated = datenum(yrs, ones(size(mons)), ones(size(days)), zeros(size(hrs)), zeros(size(mins)), zeros(size(secs)));
elseif strcmp(time_unit_string,'month')   
    datenum_truncated = datenum(yrs, mons, ones(size(days)), zeros(size(hrs)), zeros(size(mins)), zeros(size(secs)));
elseif strcmp(time_unit_string,'day')   
    datenum_truncated = datenum(yrs, mons, days, zeros(size(hrs)), zeros(size(mins)), zeros(size(secs)));
elseif strcmp(time_unit_string,'hour')
    datenum_truncated = datenum(yrs, mons, days, hrs, zeros(size(mins)), zeros(size(secs)));
elseif strcmp(time_unit_string,'minute')
    datenum_truncated = datenum(yrs, mons, days, hrs, mins, zeros(size(secs)));
else
    disp('Error: Incorrect input format or type.')
end