function VSF2D_raw = VSF2D_Calc(u,v,lat,lon,oceantime,par)
%% Info:

% USE:

%   To calculate the angle-dependent velocity structure function of any order that is of the form
%   [delta(u)^n]

% INPUT: nd x nt = matrix of dimensions number drifters x number of time steps
%       lat = nd X nt, degrees ( N = +, S = - ),
%       lon = nd x nt, degress ( E = +, W = - ), 
%       u = nd x nt, x component of velocity
%       v = nd x nt, y component of velocity
%       oceantime = row timevector (1 x nt)
%       par = structure with set of parameters 
%           par.order = order of structure function
%           par.spacesample = percent of nd to subsample
%           par.timesample = percent of nt to subsample
%           par.resamp = percent of subsampled data to pad
%           par.padding = number of grid points to pad for each resamp pt
%           par.direction = 
%               'xy' -> zonal.meridional
%               'rt' -> r, theta (coming soon)

% OUTPUT: 

%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

order = par.order;
spacesample = par.spacesample;
timesample = par.timesample;
percentresamp = par.percentresamp;
padding = par.padding;
direction = par.direction;


% subsample time
[ ~,colindx] = datasample(oceantime,timesample,2,'Replace',false); % datasample(dataset,number of samples, dimension [1=row,2=col],...)
colindx =sort(colindx);

lon = lon(:,colindx);               
lat = lat(:,colindx);
u = u(:,colindx);
v = v(:,colindx);
oceantime = oceantime(colindx);

clear colindx

% Randomly subsample

[ ~,rowindx] = datasample(lon,spacesample,1,'Replace',false); % datasample(dataset,number of samples, dimension [1=row,2=col],...)


resamp = round(percentresamp/100*length(rowindx),2);
ind = datasample(rowindx,resamp,'Replace',false);

% add in grid points plus or minus 4 grid points away from those
% originally sampled

for aa=1:padding
    indp = ind + aa.*ones(1,resamp);
    indn = ind - aa.*ones(1,resamp);
    rowindx = [rowindx indp indn];
end


clear ind indp indn resamp

rowindx =sort(rowindx);

lon = lon(rowindx,:);               
lat = lat(rowindx,:);
u = u(rowindx,:,:);
v = v(rowindx,:,:);

clear rowindx

[nd,nt] = size(lon);

%% Structure Function Values%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Define size of matrices
combos =nchoosek(nd,2);

longvel= zeros(1,combos,nt); % to add back in time add in ,nt to the third dimension slot
tranvel= zeros(1,combos,nt);
dxmap = zeros(1,combos,nt);
dymap = zeros(1,combos,nt);
dthetamap = zeros(1,combos,nt);
rd= zeros(1,combos,nt);

clear combos

if strcmp(direction,'xy') == 1

    beg = 1;
    fin = nd-1;
    for ii = 1:nd-1  
      for jj = ii+1:nd  % to put back time dependence just add in semicolons to every row
        [dx,~,~] =  m_idist(lon(ii,:),lat(ii,:),lon(jj,:),lat(ii,:));         % compute the distance between longitudes
        [dy,~,~] =  m_idist(lon(ii,:),lat(ii,:),lon(ii,:),lat(jj,:));         % compute the distance between latitudes
        [dr,~,~] =  m_idist(lon(ii,:),lat(ii,:),lon(jj,:),lat(jj,:));         % compute radial distance from one point to the other

        dx = sign(lon(jj,:) - lon(ii,:)).*dx;     % Get the correct sign for dx
        dy = sign(lat(jj,:) - lat(ii,:)).*dy;     % Get the correct sign for dy
        du = u(jj,:)-u(ii,:);           % Find the zonal velocity differences
        dv = v(jj,:)-v(ii,:);           % Find the meridional velocity differences

        r(ii,jj,:) = dr;

       % Longitudinal(sfvl) and Transverse(sfvt) sf values: 
         sfvl = ( du.*dx + dv.*dy)./dr;             % Projection of dvel onto r 
         sfvt = (-du.*dy + dv.*dx)./dr;             % Projection of dvel onto r perp

         sflvel(ii,jj,:) = sfvl.^order;                 % Adds a t vector for the lvel and tvel struc. func. values for drifters i and all others after i for all times
         sftvel(ii,jj,:) = sfvt.^order;                 % which creates an upper right triangle 2D array for each time step
         dxtemp(ii,jj,:) = dx;
         dytemp(ii,jj,:) = dy;
    %      dthetatemp(ii,jj,:) = atan2(dy,dx);
         clear sfvl sfvt dx dy du dv dr
      end

      longvel(1,beg:fin,:) = sflvel(ii,ii+1:end,:);
      tranvel(1,beg:fin,:) = sftvel(ii,ii+1:end,:);
      dxmap(1,beg:fin,:) = dxtemp(ii,ii+1:end,:);
      dymap(1,beg:fin,:) = dytemp(ii,ii+1:end,:);
    %   dthetamap(1,beg:fin,:) = dthetatemp(ii,ii+1:end,:);
      rd(1,beg:fin,:) = r(ii,ii+1:end,:);

      beg = fin+1;
      fin = fin + length(ii+2:nd);

        prev=round((ii-1)/(nd-1),2);
        if round(ii/(nd-1),2)~=prev
            disp(['PROGRESS: ', int2str(round(ii/(nd-1),2)*100), '%'])
        end
    end


    l = squeeze(longvel);     % Stores and compresses long vel data to give dimensions of velocity differences x time 
    t = squeeze(tranvel);     % Stores and compresses tran vel data to give dimensions of velocity differences x time 
    dx = squeeze(dxmap)/1000;
    dy = squeeze(dymap)/1000;
    dtheta = squeeze(dthetamap)/1000;
    rd = squeeze(rd)/1000;          % This line gives dimensions of separation distances x time, and additionally puts the distances in terms of km instead of meters

    VSF2D_raw.l = l;      % Save everything to a structure
    VSF2D_raw.t = t;    
    VSF2D_raw.dx = dx;
    VSF2D_raw.dy = dy;
    VSF2D_raw.dr = rd;
else
end


end
