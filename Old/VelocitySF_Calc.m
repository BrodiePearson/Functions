
function VSF_unbinned = VelocitySF_Calc(u,v,y,x,oceantime,par)
%% Info:

% USE:

%   To calculate the isotropic, homogenous, stationary velocity structure function of 
%   any order that is of the form [delta(u)^n] of a 2d velocity field. The
%   order is taken into account during the binning process, here only
%   delta(u) is found. If par.inhomogeneity or par.anisotropy are 'True'
%   the inhomogenous or anisotropic can be obtained.

% INPUT: np x nt = matrix of dimensions number positions x number of time steps
%       lat = np X nt, degrees ( N = +, S = - ),
%       lon = np x nt, degress ( E = +, W = - ), 
%       u = np x nt, zonal component of velocity (m/s)
%       v = np x nt, meridional component of velocity (m/s)
%       oceantime = 1 x nt row timevector (Datenumber)
%       par = structure with set of parameters for subsampling and additional features for isotropy/homogeneity 
%           par.subsample = 'True' will randomly subsample ('True'/'False')
%           par.anisotropy = 'True' keeps dx and dy instead of just r for 2D SF plotting ('True'/'False')
%           par.inhomogeneity = 'True' keeps the central location of each structure function calculation ('True'/'False')


% OPTIONAL INPUT:
%           These only needed if subsample == 'True'
%               par.percentspacesample = percent of np to subsample (integer)
%               par.percenttimesample = percent of nt to subsample (integer)
%               par.sampagain = 'True' means subsample again to get small scales ('True'/'False')
%               par.percentresamp = percent of subsampled data to pad (integer)
%               par.padding = number of grid points to pad for each resamp pt (integer)

% OUTPUT: nd x nt = matrix with dimensions of number of velocity differences x number of time steps
%       l = nd x nt longitudinal velocity differences (m/s)
%       t = nd x nt transverse velocity differences (m/s)
%       r = nd x nt radial separation distances (km)
%       oceantime = 1 x nt row timevector (Datenumber)
%
% OPTIONAL OUTPUT:
%       If par.anisotropy = 'True'
%           dx,dy = nd x nt zonal and meridional separation differences (km)
%       If par.inhomogeneity = 'True'
%           latc,lonc = nd x nt central locations of structure function calculation degrees ( N = +, S = - ) & degress ( E = +, W = - )
%       If par.subsample = 'True'
%           colindx,rowindx = subsampling indices from VelocitySF_calc


%-------------------------------------------------------------------------
% Setup
%-------------------------------------------------------------------------

[np, nt] = size(u);

%-------------------------------------------------------------------------
% Subsampling
%-------------------------------------------------------------------------

if strcmp(par.subsample,'True')
    spacesample = ceil(par.percentspacesample/100* np); 
    timesample = ceil(par.percenttimesample/100* nt);   
    
    % randomly subsample time (%nt)
    [ ~,colindx] = datasample(oceantime,timesample,2,'Replace',false); % datasample(dataset,number of samples, dimension [1=row,2=col],...)
    colindx =sort(colindx);

    x = x(:,colindx);               
    y = y(:,colindx);
    u = u(:,colindx);
    v = v(:,colindx);
    oceantime = oceantime(colindx);

    % randomly subsample space (%np)
    [ ~,rowindx] = datasample(x,spacesample,1,'Replace',false); % datasample(dataset,number of samples, dimension [1=row,2=col],...)

    if strcmp(par.sampagain,'True')
        
        % resample points to add in those nearby
        resamp = round(par.percentresamp/100*length(rowindx),2);
        ind = datasample(rowindx,resamp,'Replace',false);

        % add in grid points plus or minus 'padding' number of grid points away from those
        % resampled in the previous step

        for aa=1:par.padding
            indp = ind + aa.*ones(1,resamp);
            indn = ind - aa.*ones(1,resamp);
            rowindx = [rowindx indp indn];
        end

    end
    
    clear ind indp indn resamp

    rowindx =sort(rowindx);

    x = x(rowindx,:);               
    y = y(rowindx,:);
    u = u(rowindx,:,:);
    v = v(rowindx,:,:);

    np = size(u,1);
end


%% SF Calculation 

%*************************************************************************
%  Written by Helga Huntley
%*************************************************************************

%-------------------------------------------------------------------------
% Set up variables
%-------------------------------------------------------------------------

rtmp  = cell(1,np);   % sep distance
dxtmp = cell(1,np);   % sep in x
dytmp = cell(1,np);   % sep in y
dutmp = cell(1,np);   % zonal velocity increment
dvtmp = cell(1,np);   % meridional velocity increment

%-------------------------------------------------------------------------
% Compute distances and velocity differences
%-------------------------------------------------------------------------

parfor j = 1:np-1
    nn  = np - j;

    rtmp{j}  = m_idist(repmat(x(j,:),nn,1),...
        repmat(y(j,:),nn,1),x(j+1:np,:),y(j+1:np,:));
    
    dxtmp{j} = m_idist(repmat(x(j,:),nn,1),...
        repmat(y(j,:),nn,1),x(j+1:np,:),repmat(y(j,:),nn,1));
    dxtmp{j}(repmat(x(j,:),nn,1) < x(j+1:np,:)) = ...
        -dxtmp{j}(repmat(x(j,:),nn,1) < x(j+1:np,:));
    
    dytmp{j} = m_idist(repmat(x(j,:),nn,1),...
        repmat(y(j,:),nn,1),repmat(x(j,:),nn,1),y(j+1:np,:));
    dytmp{j}(repmat(y(j,:),nn,1) < y(j+1:np,:)) = ...
        -dytmp{j}(repmat(y(j,:),nn,1) < y(j+1:np,:));
    
    dutmp{j} = repmat(u(j,:),nn,1) - u(j+1:np,:);
    dvtmp{j} = repmat(v(j,:),nn,1) - v(j+1:np,:);
    
    prev=round((j-1)/np,2);
    if round(j/np,2)~=prev
        disp(['PROGRESS: ', int2str(round(j/np,2)*100), '%'])
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

VSF_unbinned.l = ( du.*dx + dv.*dy) ./ r;
VSF_unbinned.t = (-du.*dy + dv.*dx) ./ r;
VSF_unbinned.r =r;
VSF_unbinned.oceantime =oceantime;
VSF_unbinned.attributes = {['date modified: ' datestr(date)]};

if strcmp(par.subsample,'True')
    VSF_unbinned.colindx=colindx;
    VSF_unbinned.rowindx=rowindx;
elseif strcmp(par.inhomogeneity,'True')
    VSF_unbinned.latc= latc;
    VSF_unbinned.lonc= lonc;
elseif strcmp(par.anisotropy,'True')
    VSF_unbinned.dx = dx;
    VSF_unbinned.dy = dy;
end
    
end
