function [ul,ut,r] = rel_vel_short(x,y,u,v)

%*************************************************************************
% [ul,ut,r] = rel_vel_short(x,y,u,v)
%*************************************************************************
%  Function to retrieve pairwise relative velocities from a set of
%  positions and velocities, for use in the calculation of velocity
%  structure functions
%*************************************************************************
%
%	x,y         : positions
%	u,v         : velocities; if positions are in degrees lon and lat, 
%                 units of u and v should be m/s
% 
%   ul,ut       : longitudinal and transverse relative velocity components
%   r           : separation distance for each pair (in m, if x,y are in
%                 lon/lat)
%
%*************************************************************************
%  Written by Helga Huntley
%*************************************************************************

ntraj = size(x,1);

%-------------------------------------------------------------------------
% Set up variables
%-------------------------------------------------------------------------

rtmp  = cell(1,ntraj);   % sep distance
dxtmp = cell(1,ntraj);   % sep in x
dytmp = cell(1,ntraj);   % sep in y
dutmp = cell(1,ntraj);   % zonal velocity increment
dvtmp = cell(1,ntraj);   % meridional velocity increment

%-------------------------------------------------------------------------
% Compute distances and velocity differences
%-------------------------------------------------------------------------

parfor j = 1:ntraj-1
    nn  = ntraj - j;

    rtmp{j}  = m_idist(repmat(x(j,:),nn,1),...
        repmat(y(j,:),nn,1),x(j+1:ntraj,:),y(j+1:ntraj,:));
    
    dxtmp{j} = m_idist(repmat(x(j,:),nn,1),...
        repmat(y(j,:),nn,1),x(j+1:ntraj,:),repmat(y(j,:),nn,1));
    dxtmp{j}(repmat(x(j,:),nn,1) < x(j+1:ntraj,:)) = ...
        -dxtmp{j}(repmat(x(j,:),nn,1) < x(j+1:ntraj,:));
    
    dytmp{j} = m_idist(repmat(x(j,:),nn,1),...
        repmat(y(j,:),nn,1),repmat(x(j,:),nn,1),y(j+1:ntraj,:));
    dytmp{j}(repmat(y(j,:),nn,1) < y(j+1:ntraj,:)) = ...
        -dytmp{j}(repmat(y(j,:),nn,1) < y(j+1:ntraj,:));
    
    dutmp{j} = repmat(u(j,:),nn,1) - u(j+1:ntraj,:);
    dvtmp{j} = repmat(v(j,:),nn,1) - v(j+1:ntraj,:);
    
    prev=round((j-1)/ntraj,2);
    if round(j/ntraj,2)~=prev
        disp(['PROGRESS: ', int2str(round(j/ntraj,2)*100), '%'])
    end
end

r  = cell2mat(cat(2,rtmp(:)));
dx = cell2mat(cat(2,dxtmp(:)));
dy = cell2mat(cat(2,dytmp(:)));
du = cell2mat(cat(2,dutmp(:)));
dv = cell2mat(cat(2,dvtmp(:)));

%-------------------------------------------------------------------------
% Convert to longitudinal and transverse components
%-------------------------------------------------------------------------

ul = ( du.*dx + dv.*dy) ./ r;
ut = (-du.*dy + dv.*dx) ./ r;
