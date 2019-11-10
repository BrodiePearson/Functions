function struc = lvel_struc_func_calc(u,v,x,y,oceantime,par)
%*************************************************************************
%  Function calculate the Langrangian structure function. Output is a
%  single layer structure.
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

% make sure oceantime is a row vector
if size(oceantime,1)>size(oceantime,2)
    oceantime = oceantime';
end
% create new variables
t = repmat(oceantime,size(u,1),1);
t = datetime(t, 'ConvertFrom', 'datenum');
% switch dimensions
u = u';
v = v';
x = x';
y = y';
t = t';
nl = size(u,1);
%-------------------------------------------------------------------------
% Set up variables
%-------------------------------------------------------------------------

ttmp  = cell(1,nl);   % sep time
dutmp = cell(1,nl);   % zonal velocity increment
dvtmp = cell(1,nl);   % meridional velocity increment

%-------------------------------------------------------------------------
% Compute distances and velocity differences
%-------------------------------------------------------------------------

disp('Calculating relative velocities...')
parfor j = 1:nl-1
    nn  = nl - j;

    ttmp{j}  = days(repmat(t(j,:),nn,1) - t(j+1:nl,:));
    
    dutmp{j} = repmat(u(j,:),nn,1) - u(j+1:nl,:);

    dvtmp{j} = repmat(v(j,:),nn,1) - v(j+1:nl,:);
    
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

    prev=round((j-1)/(nl-1),2);
    if round(j/(nl-1),2)~=prev
        disp(['PROGRESS: ', int2str(round(j/(nl-1),2)*100), '%'])
    end
end

tau  = abs(cell2mat(cat(2,ttmp(:))));
du = cell2mat(cat(2,dutmp(:)));
dv = cell2mat(cat(2,dvtmp(:)));
r  = cell2mat(cat(2,rtmp(:)));
dx = cell2mat(cat(2,dxtmp(:)));
dy = cell2mat(cat(2,dytmp(:)));

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
struc.du = du;
struc.dv = dv;
struc.dul = dul;
struc.dut = dut;
struc.tau = tau;
struc.attributes = {['date modified: ' datestr(date)]};

end
