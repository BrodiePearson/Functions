function [u,v,latav,lonav,keytime] =getvelfromlatlon(lat,lon,oceantime)
%calculates the centered difference velocities in m/s for input of lat lon
%arrays with dimensions nd x nt where nd = numlocs and nt = numtimes

% Inputs:
%   lat = nd x nt array (degrees)
%   lon = nd x nt array (degrees)
%   
    for tt = 2:nt-1
        dx = m_idist(lon(:,tt-1),lat(:,tt-1),lon(:,tt+1),lat(:,tt-1));
        dx = sign(lon(:,tt+1) - lon(:,tt-1)).*dx;
        dy = m_idist(lon(:,tt-1),lat(:,tt-1),lon(:,tt-1),lat(:,tt+1));
        dy = sign(lat(:,tt+1) - lat(:,tt-1)).*dy;
        dt = 86400;
        
        u(:,tt-1) = dx/dt; % m/s
        v(:,tt-1) = dy/dt; % m/s
    end
    
    % average time, lat, and lon values
    for dd = 1:nd
        lattemp = lat(dd,:);
        lattemp = (lattemp+circshift(lattemp,2))/2;
        lattemp= lattemp(2:end-1);
        latav(dd,:) = lattemp;
        
        lontemp = lon(dd,:);
        lontemp = (lontemp+circshift(lontemp,2))/2;
        lontemp= lontemp(2:end-1);
        lonav(dd,:) = lontemp;
    end
    
    keytime = (keytime+circshift(keytime,2))/2;
    keytime = keytime(2:end-1);
    
end