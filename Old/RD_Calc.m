function RD_unbinned = RD_Calc(lat,lon,oceantime,par)
%% Info:

% USE:

%   To calculate the isotropic, homogenous, stationary erlative dispersion of a 2d velocity field

% INPUT: nd x nt = matrix of dimensions number drifters x number of time steps
%       lat = nd X nt, degrees ( N = +, S = - ),
%       lon = nd x nt, degress ( E = +, W = - ), 
%       u = nd x nt, x component of velocity
%       v = nd x nt, y component of velocity
%       oceantime = 1 x nt row timevector 
%       par = structure with set of parameters 
%           par.percentspacesample = percent of nd to subsample (ex: 10,50)
%           par.percenttimesample = percent of nt to subsample
%           par.resamp = percent of subsampled data to pad
%           par.padding = number of grid points to pad for each resamp pt


% OUTPUT: nd x nt = matrix with dimensions of number of velocity
% differences by number of time steps
%       l = nd x nt longitudinal structure function values (m/s)^2
%       t = nd x nt transverse structure function values (m/s)^2
%       r = nd x nt separation distances (km)
%       oceantime = 1 x nt row timevector

%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nd,nt] = size(lat);

%nd = 10;

counter = 1; % variable to keep track of the progress of the code

% Define size of matrices
combos =nchoosek(nd,2);

zd = zeros(combos,nt);
zm = zeros(combos,nt);
s = zeros(combos,nt);

    parfor tt = 1:nt % calculate separately for each time step in parallel

        %Assign space for variables
        lonvel= zeros(1,combos); 
        tranvel= zeros(1,combos);
        rd= zeros(1,combos);
        r = zeros(nd-1,nd);
        sflvel = zeros(nd-1,nd);
        sftvel = zeros(nd-1,nd);

        beg = 1;
        fin = nd-1;
        for ii = 1:nd-1  
          for jj = ii+1:nd  
            [dx,~,~] =  m_idist(lon(ii,tt),lat(ii,tt),lon(jj,tt),lat(ii,tt));         % compute the distance between longitudes
            [dy,~,~] =  m_idist(lon(ii,tt),lat(ii,tt),lon(ii,tt),lat(jj,tt));         % compute the distance between latitudes

            r(ii,jj) =  m_idist(lon(ii,tt),lat(ii,tt),lon(jj,tt),lat(jj,tt));         % compute radial distance from one point to the other
            sflvel(ii,jj) = dx.^2;                 % Adds a t vector for the lvel and tvel struc. func. values for drifters i and all others after i for all times
            sftvel(ii,jj) = dy.^2;                 % which creates an upper right triangle 2D array for each time step

          end

          lonvel(1,beg:fin) = sflvel(ii,ii+1:end);
          tranvel(1,beg:fin) = sftvel(ii,ii+1:end);
          rd(1,beg:fin) = r(ii,ii+1:end);

          beg = fin+1;
          fin = fin + length(ii+2:nd);

        end


        zd(:,tt) = squeeze(lonvel);     % Stores and compresses long vel data to give dimensions of velocity differences x time 
        zm(:,tt) = squeeze(tranvel);     % Stores and compresses tran vel data to give dimensions of velocity differences x time 

        s(:,tt) = squeeze(rd)/1000;          % This line gives dimensions of separation distances x time, and additionally puts the distances in terms of km instead of meters



        prev=round((tt-1)/nt,2);
        if round(tt/nt,2)~=prev
            disp(['PROGRESS Calculating: ', int2str(round(tt/nt,2)*100), '%'])
        end

    end

    RD_unbinned.oceantime =oceantime;
    RD_unbinned.zd = zd;
    RD_unbinned.md = md;
    RD_unbinned.r= r;
    
    RD_unbinned.DM = ['date modified: ' datestr(date)];
    RD_unbinned.attributes = {['space subsample: ' num2str(spacesample)];['time subsample: ' num2str(timesample)]; ['percent resampled: ' num2str(percentresamp) '%'];[ 'padding: ' num2str(padding)]};


end