function dvel = deltavel(u,v,lat,lon, oceantime)
%% Info:


%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nd, nt] = size(u);

%% SF Calculation 
   
    counter = 1; % variable to keep track of the progress of the code
    
    % Define size of matrices
    combos =nchoosek(nd,2);
    
    l = zeros(combos,nt);
    t = zeros(combos,nt);
    s = zeros(combos,nt);
    
    
    parfor tt = 1:nt % calculate separately for each time step

    
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
      for jj = ii+1:nd  % to put back time dependence just add in semicolons to every row
        [dx,~,~] =  m_idist(lon(ii,tt),lat(ii,tt),lon(jj,tt),lat(ii,tt));         % compute the distance between longitudes
        [dy,~,~] =  m_idist(lon(ii,tt),lat(ii,tt),lon(ii,tt),lat(jj,tt));         % compute the distance between latitudes
        [dr,~,~] =  m_idist(lon(ii,tt),lat(ii,tt),lon(jj,tt),lat(jj,tt));         % compute radial distance from one point to the other

        dx = sign(lon(jj,tt) - lon(ii,tt)).*dx;     % Get the correct sign for dx
        dy = sign(lat(jj,tt) - lat(ii,tt)).*dy;     % Get the correct sign for dy

        du = u(jj,tt)-u(ii,tt);           % Find the zonal velocity differences
        dv = v(jj,tt)-v(ii,tt);           % Find the meridional velocity differences

        r(ii,jj) = dr;

       % Longitudinal(sfvl) and Transverse(sfvt) sf values: 
         sfvl = ( du.*dx + dv.*dy)./dr;             % Projection of dvel onto r 
         sfvt = (-du.*dy + dv.*dx)./dr;             % Projection of dvel onto r perp

         sflvel(ii,jj) = sfvl;                 % Adds a t vector for the lvel and tvel struc. func. values for drifters i and all others after i for all times
         sftvel(ii,jj) = sfvt;                 % which creates an upper right triangle 2D array for each time step

         
         counter = counter+1; % keep track of progress
      end

      lonvel(1,beg:fin) = sflvel(ii,ii+1:end);
      tranvel(1,beg:fin) = sftvel(ii,ii+1:end);
      rd(1,beg:fin) = r(ii,ii+1:end);

      beg = fin+1;
      fin = fin + length(ii+2:nd);

    end


    l(:,tt) = squeeze(lonvel);     % Stores and compresses long vel data to give dimensions of velocity differences x time 
    t(:,tt) = squeeze(tranvel);     % Stores and compresses tran vel data to give dimensions of velocity differences x time 

    s(:,tt) = squeeze(rd)/1000;          % This line gives dimensions of separation distances x time, and additionally puts the distances in terms of km instead of meters

    
    prev=round((tt-1)/nt,2);
    if round(tt/nt,2)~=prev
        disp(['PROGRESS Calculating: ', int2str(round(tt/nt,2)*100), '%'])
    end
    
    end
    
    
    dvel.DM = ['date modified: ' datestr(date)];
    dvel.dl = l;
    dvel.dt = t;
    dvel.dr = s;
    dvel.oceantime = oceantime;


end