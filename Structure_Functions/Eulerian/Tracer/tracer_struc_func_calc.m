function struc = tracer_struc_func_calc(beta,lon,lat,oceantime,par)

%*************************************************************************
%  struc = tracer_struc_func_calc(x,y,u,v,oceantime,par)
%*************************************************************************
%  Function to retrieve pairwise tracer differences from a set of
%  positions and tracer values, for use in the calculation of tracer and
%  bblended structure functions. Output is single layer structure.
%
% nl = no. locations, nt = no. times, np = no. pairs
%*************************************************************************
%   Required Input
%-------------------------------------------------------------------------
%	x,y         : positions [nl,nt] (lon,lat in deg. or m)
%	theta       : tracer [nl,nt] ([theta])
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
%   dbeta       : tracer incremements [np,nt] ([theta])
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
if nargin < 5
    par.empty = 'True';
end

nl = size(lon,1);

%-------------------------------------------------------------------------
% Set up variables
%-------------------------------------------------------------------------

rtmp  = cell(1,nl);   % sep distance
dxtmp = cell(1,nl);   % sep in x
dytmp = cell(1,nl);   % sep in y
dbetatmp = cell(1,nl);   % tracer increment

%-------------------------------------------------------------------------
% Compute distances and velocity differences
%-------------------------------------------------------------------------

disp('Calculating tracer increments.')
parfor j = 1:nl-1
    nn  = nl - j

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
    
    dbetatmp{j} = repmat(beta(j,:),nn,1) - beta(j+1:nl,:);

    prev=round((j-1)/(nl-1),2);
    if round(j/(nl-1),2)~=prev
        disp(['PROGRESS: ', int2str(round(j/(nl-1),2)*100), '%'])
    end
end

r  = cell2mat(cat(2,rtmp(:)));
dx = cell2mat(cat(2,dxtmp(:)));
dy = cell2mat(cat(2,dytmp(:)));
dbeta = cell2mat(cat(2,dbetatmp(:)));

%-------------------------------------------------------------------------
% Save data to structure
%-------------------------------------------------------------------------
disp('Saving data to structure.')

struc.dx = dx./1000;
struc.dy = dy./1000;
struc.dbeta = dbeta;
struc.r = r./1000;
struc.oceantime =oceantime;
struc.par = par;
end
