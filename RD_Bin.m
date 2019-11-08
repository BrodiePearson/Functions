function VSF_binned = RD_Bin(zd,md,r,oceantime,par)
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




if strcmp(par.direction,'xy') == 1
    %% Setup

    [nd,nt] = size(x);

    %% Structure Function Binning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % create a logarithmically spaced bin vector that increases by a power of 10 every four points
    % Bins(1) = min(r(:));
    Bins(1) = .5;
    kk = 2;
    while(Bins(kk-1) <= ceil(max(r(:))))
        Bins(kk) = Bins(kk-1)*10^((kk-2)*par.n);  
        kk = kk+1;
    end


    % Time-stamped structure functions
    for tt = 1:nt

        for ii = 1:length(Bins)-1

            ind = find(r(:,tt) <= Bins(ii+1) & r(:,tt) > Bins(ii));   % Isolate the indices of data that fall into a given bin

            longvel = zd(:,tt);
            tranvel = md(:,tt);

            sflv = longvel(ind);                            % Find long. sf values associated with the indices from ^
            sftv = tranvel(ind);                            % Find tran. sf values associated with the indices from ^
            tsf(ii) = mean(sftv(~isnan(sftv)));                        % Compute the average of all data in the tran. sf bin and save to vector for plotting
            lsf(ii) = mean(sflv(~isnan(sflv)));                        % Compute the average of all data in the long. sf bin and save to vector for plotting
            counttsf(ii) = length(sftv(~isnan(sftv)));
            countlsf(ii) = length(sflv(~isnan(sflv)));


            prev=round((ii-1)/(length(Bins)-1),2);
            if round(ii/(length(Bins)-1),2)~=prev
                disp(['PROGRESS Binning: ', int2str(round(ii/(length(Bins)-1),2)*100), '%'])
            end

        end

         % Add info to matrix

         TSF(:,tt) = tsf;
         LSF(:,tt) = lsf;
         COUNTTSF(:,tt) = counttsf;
         COUNTLSF(:,tt) = countlsf;

    end

    % Take Ensemble Mean __> sum(....,2) sums along the rows to give a column
    % vector of the row sums

    sumtsf = nansum(TSF.*COUNTTSF,2);
    tottsf = nansum(COUNTTSF,2);
    enstsf = sumtsf./tottsf;

    sumlsf = nansum(LSF.*COUNTLSF,2);
    totlsf = nansum(COUNTLSF,2);
    enslsf = sumlsf./totlsf;

    % Center Bins
    disp('Centering Bins...')
    rmid(1)=par.minbin;
    for jj = 2:length(Bins)-1
        rmid(jj) = (Bins(jj)+Bins(jj+1))/2;     % Declare the average distance between bins as the representative displacement for a given bin ( centering )
    end
  
    
    order = par.order;
    spacesample = ceil(par.percentspacesample/100* nd);
    timesample = ceil(par.percenttimesample/100* nt);
    percentresamp = par.percentresamp;
    padding = par.padding;
    direction = par.direction;

    
    VSF_binned.enszd = enslsf;
    VSF_binned.ensmd = enstsf;
    VSF_binned.r = rmid;
    VSF_binned.DM = ['date modified: ' datestr(date)];
    VSF_binned.attributes = {['space subsample: ' num2str(spacesample)];['time subsample: ' num2str(timesample)]; ['percent resampled: ' num2str(percentresamp) '%'];[ 'padding: ' num2str(padding)]; ['binning power: ' num2str(par.n)]};
    
end


end