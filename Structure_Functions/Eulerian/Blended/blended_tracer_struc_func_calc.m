function struc = blended_tracer_struc_func_calc(x,y,u,v,beta,oceantime,par)
%*************************************************************************
%  Function to retrieve pairwise relative velocities from a set of
%  positions and velocities, for use in the calculation of velocity
%  structure functions. Output is single layer structure.
%
% nl = no. locations, nt = no. times, np = no. pairs
%*************************************************************************
%  Syntax:
%-------------------------------------------------------------------------
%   struc = blended_tracer_struc_func_calc(x,y,u,v,beta,oceantime,Property,Value,...)
%*************************************************************************
%   Required Input
%-------------------------------------------------------------------------
%	x,y         : positions [nl,nt] (lon,lat in deg. or m)
%	u,v         : meridinoal and zonal velocity [nl,nt] (m/s)
%   oceantime   : time [1,nt] (datenumber)
%-------------------------------------------------------------------------
%   Optional Input
%------------------------------------------------------------------------- 
%   anisotropy  : if set to 'True' saves and returns dx,dy
%   homogeneity : if set to 'True' saves and returns latc, lonc
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
%  Written by Jenna Pearson
%*************************************************************************

%-------------------------------------------------------------------------
% Error Checks and Parse Optional Inputs
%-------------------------------------------------------------------------
% narginchk(6,inf)
% % assert(isequal(size(x),size(y),size(u),size(v),size(beta))==1,'Error: Dimensions of x, y, u, v, and beta must agree. Check inputs, or use [x,y] = meshgrid(x,y) if necessary.')
% % assert(isequal(size(x,2),length(oceantime)),'Error: Length of variable oceantime must match second dimension of x, y, u, v, and beta.') 
% 
% tmp = strncmpi(varargin,'inhomogeneity',3); 
% if any(tmp)
%    inhomogeneity = varargin{find(tmp)+1}; 
%    tmp(find(tmp)+1)=1; 
%    varargin = varargin(~tmp); 
%    assert(islogical(inhomogeneity),'Input error: Inhomogeneity option must be true or false. Default is true.') 
% end
% 
% tmp = strncmpi(varargin,'anisotropy',3); 
% if any(tmp)
%    anisotropy = varargin{find(tmp)+1}; 
%    tmp(find(tmp)+1)=1; 
%    varargin = varargin(~tmp); 
%    assert(islogical(anisotropy),'Input error: Anisotropy option must be true or false. Default is true.') 
% end

%-------------------------------------------------------------------------
% Set up variables
%-------------------------------------------------------------------------
nl = size(x,1);
rtmp  = cell(1,nl);   % sep distance
dxtmp = cell(1,nl);   % sep in x
dytmp = cell(1,nl);   % sep in y
dutmp = cell(1,nl);   % zonal velocity increment
dvtmp = cell(1,nl);   % meridional velocity increment
dbetatmp = cell(1,nl);   % tracer increment
% 
% if inhomogeneity
%     latctmp = cell(1,nl); % mean lat of sf calcs
%     lonctmp = cell(1,nl); % mean lon of sf calcs
% end


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
    dbetatmp{j} = repmat(beta(j,:),nn,1) - beta(j+1:nl,:);

end

r  = cell2mat(cat(2,rtmp(:)));
dx = cell2mat(cat(2,dxtmp(:)));
dy = cell2mat(cat(2,dytmp(:)));
du = cell2mat(cat(2,dutmp(:)));
dv = cell2mat(cat(2,dvtmp(:)));
dbeta = cell2mat(cat(2,dbetatmp(:)));

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
struc.dx = dx./1000;
struc.dy = dy./1000;
struc.dbeta=dbeta;
struc.r = r./1000;
struc.oceantime =oceantime;
end
