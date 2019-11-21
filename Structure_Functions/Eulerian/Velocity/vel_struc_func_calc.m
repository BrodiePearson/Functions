function struc = vel_struc_func_calc(x,y,u,v,oceantime,par)

%*************************************************************************
%  struc = vel_struc_func_calc(x,y,u,v,oceantime,par)
%*************************************************************************
%  Function to retrieve pairwise relative velocities from a set of
%  positions and velocities, for use in the calculation of velocity
%  structure functions. Output is single layer structure.
%
% nl = no. locations, nt = no. times, np = no. pairs
%*************************************************************************
%   Required Input
%-------------------------------------------------------------------------
%	x,y         : positions [nl,nt] (lon,lat in deg. or m)
%	u,v         : meridinoal and zonal velocity [nl,nt] (m/s)
%   oceantime   : time [1,nt] (datenumber)
%-------------------------------------------------------------------------
%   Optional Input
%-------------------------------------------------------------------------
%   par         : a structure containing the following variables
%   
%   anisotropy  : if set to 'True' saves dx,dy
%   homogeneity : if set to 'True' saves latc, lonc
%-------------------------------------------------------------------------
%   Default Output
%-------------------------------------------------------------------------
%   ul,ut       : longitudinal and transverse relative velocity components
%                 [np,nt] (m/s)
%   r           : separation distance for each pair [np,nt] (km)
%   oceantime   : time [1,nt] (datenumber)
%   attributes  : structure containing the date modified
%-------------------------------------------------------------------------
%   Optional Output
%-------------------------------------------------------------------------
%   dx,dy       : meridional and zone separations [np,nt] (m)
%   latc,lonc   : average position of sf values [np,nt] (lat lon in deg.)
%*************************************************************************
%  Written by Jenna Pearson, and adapted from Helga Huntley
%*************************************************************************

% check to make sure par is defined
if nargin < 6
    par.empty = 'True';
end

nl = size(x,1);

%-------------------------------------------------------------------------
% Set up variables
%-------------------------------------------------------------------------

rtmp  = cell(1,nl);   % sep distance
dxtmp = cell(1,nl);   % sep in x
dytmp = cell(1,nl);   % sep in y
dutmp = cell(1,nl);   % zonal velocity increment
dvtmp = cell(1,nl);   % meridional velocity increment

if  isfield(par,'inhomogeneity') && strcmp(par.inhomogeneity,'True') 
    latctmp = cell(1,nl); % mean lat of sf calcs
    latctmp = cell(1,nl); % mean lat of sf calcs
end

%-------------------------------------------------------------------------
% Compute distances and velocity differences
%-------------------------------------------------------------------------

disp('Calculating relative velocities...')
parfor j = 1:nl-1
    nn  = nl - j;

    rtmp{j}  = m_idist(repmat(x(j,:),nn,1),...
        repmat(y(j,:),nn,1),x(j+1:nl,:),y(j+1:nl,:));
    
    dxtmp{j} = m_idist(repmat(x(j,:),nn,1),...
        repmat(y(j,:),nn,1),x(j+1:nl,:),repmat(y(j,:),nn,1));
    dxtmp{j}(repmat(x(j,:),nn,1) < x(j+1:nl,:)) = ...
        -dxtmp{j}(repmat(x(j,:),nn,1) < x(j+1:nl,:));
    
    dytmp{j} = m_idist(repmat(x(j,:),nn,1),...
        repmat(y(j,:),nn,1),repmat(x(j,:),nn,1),y(j+1:nl,:));
    dytmp{j}(repmat(y(j,:),nn,1) < y(j+1:nl,:)) = ...
        -dytmp{j}(repmat(y(j,:),nn,1) < y(j+1:nl,:));
    
    dutmp{j} = repmat(u(j,:),nn,1) - u(j+1:nl,:);
    dvtmp{j} = repmat(v(j,:),nn,1) - v(j+1:nl,:);
    
    if  isfield(par,'inhomogeneity') && strcmp(par.inhomogeneity,'True') 
        %finish this
    end
%     
%     if mod(round(j/nl,2)*100,25)==0
%         disp([int2str(round(j/nl,2)*100) '% Complete'])
%     end
end

r  = cell2mat(cat(2,rtmp(:)));
dx = cell2mat(cat(2,dxtmp(:)));
dy = cell2mat(cat(2,dytmp(:)));
du = cell2mat(cat(2,dutmp(:)));
dv = cell2mat(cat(2,dvtmp(:)));

%-------------------------------------------------------------------------
% Convert to longitudinal and transverse components
%-------------------------------------------------------------------------
disp('Converting to ul and ut...')
ul = ( du.*dx + dv.*dy) ./ r;
ut = (-du.*dy + dv.*dx) ./ r;

%-------------------------------------------------------------------------
% Save data to structure
%-------------------------------------------------------------------------
disp('Saving data to structure...')
struc.ul = ul;
struc.ut = ut;
struc.r = r./1000;
struc.oceantime =oceantime;
struc.attributes = {['date modified: ' datestr(date)]};
struc.par = par;

if isfield(par,'inhomogeneity') && strcmp(par.inhomogeneity,'True') % need to add this part in 
    struc.latc= latc;
    struc.lonc= lonc;
end

if isfield(par,'anisotropy') && strcmp(par.anisotropy,'True')
    struc.dx = dx;
    struc.dy = dy;
end
end
