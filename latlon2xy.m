function [x,y] = latlon2xy(lat,lon,units)

if nargin < 3
    units = 'km';
end

h = 0;
h0 = 0;
spheroid = referenceEllipsoid('GRS 80',units);


[x,y,~] = geodetic2ned(lat,lon,h,nanmean(lat,1),nanmean(lon,1),h0,spheroid);
end