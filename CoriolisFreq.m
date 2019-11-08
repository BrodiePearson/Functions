function [f] = CoriolisFreq(lat)

% Inputs:
%   lat = latitude in degrees
% Outputs:
%   f = coriolisf(lat) returns the Coriolis frequency in radians per second for point(s) at latitude(s) lat.  



assert(max(abs(lat(:)))<=90,'Latitude value(s) out of realistic bounds. Check inputs.')


f = 2*(7.2921e-5).*sind(lat); 
end