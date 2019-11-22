function struc = vel_struc_func_calc(u,v,lon,lat,oceantime,par)

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

nl = size(lon,1);

%-------------------------------------------------------------------------
% Set up variables
%-------------------------------------------------------------------------

rtmp  = cell(1,nl);   % sep distance
dxtmp = cell(1,nl);   % sep in x
dytmp = cell(1,nl);   % sep in y
dutmp = cell(1,nl);   % zonal velocity increment
dvtmp = cell(1,nl);   % meridional velocity increment


if isfield(par,'inhomogeneity') && strcmp(par.inhomogeneity,'True') 
    latctmp = cell(1,nl);
    lonctmp = cell(1,nl); 
end


%-------------------------------------------------------------------------
% Compute distances and velocity differences
%-------------------------------------------------------------------------

disp('Calculating relative velocities...')
parfor j = 1:nl-1
    nn  = nl - j;

    rtmp{j}  = m_idist(repmat(lon(j,:),nn,1),...
        repmat(lat(j,:),nn,1),lon(j+1:nl,:),lat(j+1:nl,:));
    
    dxtmp{j} = m_idist(repmat(lon(j,:),nn,1),...
        repmat(lat(j,:),nn,1),lon(j+1:nl,:),repmat(lat(j,:),nn,1));
    dxtmp{j}(repmat(lon(j,:),nn,1) < lon(j+1:nl,:)) = ...
        -dxtmp{j}(repmat(lon(j,:),nn,1) < lon(j+1:nl,:));
    
    dytmp{j} = m_idist(repmat(lon(j,:),nn,1),...
        repmat(lat(j,:),nn,1),repmat(lon(j,:),nn,1),lat(j+1:nl,:));
    dytmp{j}(repmat(lat(j,:),nn,1) < lat(j+1:nl,:)) = ...
        -dytmp{j}(repmat(lat(j,:),nn,1) < lat(j+1:nl,:));
    
    dutmp{j} = repmat(u(j,:),nn,1) - u(j+1:nl,:);
    dvtmp{j} = repmat(v(j,:),nn,1) - v(j+1:nl,:);
    

    if isfield(par,'inhomogeneity') && strcmp(par.inhomogeneity,'True') 
        latctmp{j} = (repmat(lat(j,:),nn,1) + lat(j+1:nl,:))./2;
        lonctmp{j} = (repmat(lon(j,:),nn,1) + lon(j+1:nl,:))./2;
    end

    prev = mod(round((j-1)/nl,2)*100,25);
    if mod(round(j/nl,2)*100,25)~=prev
        disp([int2str(round(j/nl,2)*100) '% Complete'])
    end
end

r  = cell2mat(cat(2,rtmp(:)));
dx = cell2mat(cat(2,dxtmp(:)));
dy = cell2mat(cat(2,dytmp(:)));
du = cell2mat(cat(2,dutmp(:)));
dv = cell2mat(cat(2,dvtmp(:)));

if isfield(par,'inhomogeneity') && strcmp(par.inhomogeneity,'True') 
    latc = cell2mat(cat(2,latctmp(:)));
    lonc = cell2mat(cat(2,lonctmp(:)));
end

%-------------------------------------------------------------------------
% Convert to longitudinal and transverse components
%-------------------------------------------------------------------------
disp('Converting to ul and ut...')
dul = ( du.*dx + dv.*dy) ./ r;
dut = (-du.*dy + dv.*dx) ./ r;

%-------------------------------------------------------------------------
% Save data to structure
%-------------------------------------------------------------------------
disp('Saving data to structure...')
struc.dul = dul;
struc.dut = dut;
struc.r = r./1000;
struc.oceantime =oceantime;
struc.attributes = {['date modified: ' datestr(date)]};
struc.par = par;

if isfield(par,'inhomogeneity') && strcmp(par.inhomogeneity,'True') 
    struc.latc= latc;
    struc.lonc= lonc;
end

if isfield(par,'anisotropy') && strcmp(par.anisotropy,'True')
    struc.dx = dx;
    struc.dy = dy;
end
end
