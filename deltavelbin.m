function dvelbin = deltavelbin(l,t,r)
%% Info:

% USE:

%   To bin and bootstrap the calculated isotropic, homogenous, stationary 
%   velocity structure function output from VSF_Calc

% INPUT: nd x nt = matrix of dimensions number drifters x number of time steps
%       l = nd x nt longitudinal structure function values (m/s)^2
%       t = nd x nt transverse structure function values (m/s)^2
%       r = nd x nt separation distances (km)
%       oceantime = 1 x nt row vector of time steps
%       par = structure with set of parameters 
%           par.n = exponent of bins(.25 --> ^ by power of 10 every 4 bins)
%           par.direction = 
%               'xy' -> zonal.meridional 
%               'lt' -> longitudinal, transverse

% OUTPUT: nd x nt = matrix with dimensions of number of velocity
% differences by number of time steps
%       l = nd x nt longitudinal structure function values (m/s)^2
%       t = nd x nt transverse structure function values (m/s)^2
%       r = nd x nt separation distances (km)
%       oceantime = 1 x nt row vector of time steps


    %% Setup

    [nd,nt] = size(l);
  
    TSF = {};
    LSF = {};

    %% Structure Function Binning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % create a logarithmically spaced bin vector that increases by a power of 10 every four points
    Bins(1) = min(r(:));
    kk = 2;
    while(Bins(kk-1) <= ceil(max(r(:))))
        Bins(kk) = Bins(kk-1)*10^((kk-1)*.3);  
        kk = kk+1;
    end
    % Time-stamped structure functions
    for tt = 1:nt

        for ii = 1:length(Bins)-1

            ind = find(r(:,tt) <= Bins(ii+1) & r(:,tt) > Bins(ii));   % Isolate the indices of data that fall into a given bin

            longvel = l(:,tt);
            tranvel = t(:,tt);

            sflv = longvel(ind);                           % Find long. sf values associated with the indices from ^
            sftv = tranvel(ind) ;                           % Find tran. sf values associated with the indices from ^
        

         % Add info to matrix
         
       
            if isempty(TSF)==1 && isempty(sftv) == 0
                TSF{ii} = sftv;
                LSF{ii} = sflv;
            elseif isempty(sftv) == 0 && ii>1     
                TSF{ii} = [TSF{ii}; sftv];
                LSF{ii} = [TSF{ii}; sflv];

            end
        end
        
        prev=round((tt-1)/(nt),2);
        if round(tt/(nt),2)~=prev
            disp(['PROGRESS Binning: ', int2str(round(tt/(nt),2)*100), '%'])
        end
    end


    % Center Bins
    disp('Centering Bins...')
    for jj = 2:length(Bins)-1
        rmid(jj) = (Bins(jj)+Bins(jj+1))/2;     % Declare the average distance between bins as the representative displacement for a given bin ( centering )
    end

    dvelbin.dl = TSF;
    dvelbin.dt = LSF;
    dvelbin.r = rmid;
end