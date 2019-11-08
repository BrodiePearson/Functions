function NSF_unbinned = NutrientSF_Calc(nut,lat,lon,oceantime,par)
%% Info:

% USE:

%   To calculate the isotropic, homogenous, stationary velocity structure function of 
%   any order that is of the form [delta(u)^n] of a 2d velocity field. The
%   order is taken into account during the binning process, here only
%   delta(u) is found.

% INPUT: nd x nt = matrix of dimensions number drifters x number of time steps
%       lat = nd X nt, degrees ( N = +, S = - ),
%       lon = nd x nt, degress ( E = +, W = - ), 
%       nut = nd x nt, x nutrient
%       oceantime = 1 x nt row timevector 
%       par = structure with set of parameters
%           par.subsample = 'True' with subsampling, 'False' without
%           These only needed if subsample == 'True'
%               par.percentspacesample = percent of nd to subsample (ex: 10,50)
%               par.percenttimesample = percent of nt to subsample
%               par.resamp = percent of subsampled data to pad
%               par.padding = number of grid points to pad for each resamp pt
%           par.direction = 
%               'xy' -> zonal.meridional 
%               'lt' -> longitudinal, transverse

% OUTPUT: nd x nt = matrix with dimensions of number of velocity
% differences by number of time steps
%       n = nd x nt nutrient structure function values,unbinned
%       r = nd x nt separation distances (km)
%       oceantime = 1 x nt row timevector

%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nd, nt] = size(nut);

%% SF Calculation 

if strcmp(par.direction,'lt') ==1

% Define size of matrices
combos =nchoosek(nd,2);
    
n = zeros(combos,nt);
s = zeros(combos,nt);
    
parfor tt = 1:length(oceantime)
    nuttemp= zeros(1,combos); 
    rd= zeros(1,combos);
    latcen= zeros(1,combos);
    loncen= zeros(1,combos);
    deltax = zeros(1,combos);
    deltay = zeros(1,combos);
    r = zeros(nd-1,nd);
    nutsf = zeros(nd-1,nd);
    latcentemp = zeros(nd-1,nd);
    loncentemp = zeros(nd-1,nd);
    dxmap = zeros(nd-1,nd);
    dymap = zeros(nd-1,nd);

    beg = 1;
    fin = nd-1;
    for ii = 1:nd-1  
      for jj = ii+1:nd  % to put back time dependence just add in semicolons to every row
        [dxtemp,~,~] =  m_idist(lon(ii,tt),lat(ii,tt),lon(jj,tt),lat(ii,tt));         % compute the distance between longitudes
        [dytemp,~,~] =  m_idist(lon(ii,tt),lat(ii,tt),lon(ii,tt),lat(jj,tt));         % compute the distance between latitudes
        [dr,~,~] =  m_idist(lon(ii,tt),lat(ii,tt),lon(jj,tt),lat(jj,tt));         % compute radial distance from one point to the other

        latcentemp(ii,jj) = mean([lat(ii,tt) lat(jj,tt)]);
        loncentemp(ii,jj) = mean([lon(ii,tt) lon(jj,tt)]);

        dxtemp = sign(lon(jj,tt) - lon(ii,tt)).*dxtemp;     % Get the correct sign for dx
        dytemp = sign(lat(jj,tt) - lat(ii,tt)).*dytemp;     % Get the correct sign for dy
        dn = nut(jj,tt)-nut(ii,tt);
        r(ii,jj) = dr;          % Projection of dvel onto r perp

        nutsf(ii,jj) = dn;  
        dxmap(ii,jj) = dxtemp;
        dymap(ii,jj) = dytemp;
      end

      nuttemp(1,beg:fin) = nutsf(ii,ii+1:end);
      rd(1,beg:fin) = r(ii,ii+1:end);
      latcen(1,beg:fin) = latcentemp(ii,ii+1:end);
      loncen(1,beg:fin) = loncentemp(ii,ii+1:end);
      deltax(1,beg:fin) = dxmap(ii,ii+1:end);
      deltay(1,beg:fin) = dymap(ii,ii+1:end);
      beg = fin+1;
      fin = fin + length(ii+2:nd);
    end

    n(:,tt) = squeeze(nuttemp);
    s(:,tt) = squeeze(rd)/1000;           
    latc(:,tt) = squeeze(latcen);
    lonc(:,tt) = squeeze(loncen);
    dx(:,tt) = squeeze(deltax)/1000;
    dy(:,tt) = squeeze(deltay)/1000;

        prev=round((tt-1)/nt,2);
        if round(tt/nt,2)~=prev
            disp(['PROGRESS: ', int2str(round(tt/nt,2)*100), '%'])
        end

end

    NSF_unbinned.n = n;
    NSF_unbinned.r= s;
    NSF_unbinned.latc= latc;
    NSF_unbinned.lonc= lonc;
    NSF_unbinned.dx = dx;
    NSF_unbinned.dy = dy;
    NSF_unbinned.oceantime =oceantime;
    NSF_unbinned.attributes = {['date modified: ' datestr(date)]};
    
    if strcmp(par.subsample,'True') == 1
        NSF_unbinned.colindx=colindx;
        NSF_unbinned.rowindx=rowindx;
    end
end
end