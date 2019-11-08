function NSF_binned = NutrientSF_bin(n,r,oceantime,par)
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


if strcmp(par.direction,'lt') == 1
    
    %% Setup

    n = n.^par.order;
    
    [nd,nt] = size(n);

    %% Structure Function Binning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % create a logarithmically spaced bin vector that increases by a power of 10 every four points
    %% Structure Function Binning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Bins(1) = 0;
    kk = 2;

    while(10^((kk-3)*.25) <= ceil(max(r(:))))
        Bins(kk-1) = 10^((kk-3)*.25);                       % Logarithmically spaced bin vector that increases by a power of 10 every four points
        kk = kk+1;
    end



    % Time-stamped structure functions
    for tt = 1:length(oceantime)

        for ii = 1:length(Bins)-1

            ind = find(r(:,tt) < Bins(ii+1) & r(:,tt) >= Bins(ii));   % Isolate the indices of data that fall into a given bin

            nuttemp = n(:,tt);
            nuttemp = nuttemp(ind);                            % Find long. sf values associated with the indices from ^
            nsf(ii) = mean(nuttemp(~isnan(nuttemp)));                        % Compute the average of all data in the long. sf bin and save to vector for plotting
            countnsf(ii) = length(nuttemp(~isnan(nuttemp)));
            NSFmed{ii,tt} = nuttemp;
        end

         % Add info to matrix
         NSF(:,tt) = nsf;
         COUNTNSF(:,tt) = countnsf;

            prev=round((tt-1)/length(oceantime),2);
            if round(tt/length(oceantime),2)~=prev
                disp(['PROGRESS: ', int2str(round(tt/length(oceantime),2)*100), '%'])
            end

    end

    % Take Ensemble Mean along the columns 

    sumnsf = nansum(NSF.*COUNTNSF,2);
    totnsf = nansum(COUNTNSF,2);
    ensnsf = sumnsf./totnsf;

    % Center Bins
    disp('Centering Bins...')
    for jj = 1:length(Bins)-1
        rmid(jj) = (Bins(jj)+Bins(jj+1))/2;     % Declare the average distance between bins as the representative displacement for a given bin ( centering )
    end

    %% Bootstrapping
    
    disp('Boostrapping...')

    [sz,~] = size(NSF);

    for kk = 1:sz
        if isnan(nanmean(NSF(kk,:))) == 0
            [btcin,] = bootci(10000,@nanmean,NSF(kk,:));    % Bootstrap LSF to find confidence intervals
            nci(kk,1) = btcin(1,1);                       % Save lower confidence interval value
            nci(kk,2) = btcin(2,1);                       % Save upper confidence interval value
        else
            nci(kk,1) = NaN;                       % Save lower confidence interval value
            nci(kk,2) = NaN; 
        end
        
        prev = round(kk-1/sz,2);
        
        if round(kk/sz,2)~=prev
            disp([int2str(round(kk/sz,2)*100) '%'])
        end
    end
    




    %% Write Data

    NSF_binned.nsf = ensnsf;
    NSF_binned.countnsf = COUNTNSF;
    NSF_binned.countnsftotals = sum(COUNTNSF,2);

    NSF_binned.NSF = NSF;
    
    NSF_binned.NSFmed = NSFmed;
    
    NSF_binned.oceantime = oceantime;

    NSF_binned.nCIL= nci(:,1);
    NSF_binned.nCIU= nci(:,2);
    
    NSF_binned.n = ensnsf;
    NSF_binned.r = rmid;
    NSF_binned.attributes = ['date modified: ' datestr(date)];
    
    if strcmp(par.subsample,'True') == 1
        NSF_binned.colindx=par.colindx;
        NSF_binned.rowindx=par.rowindx;
    end
    
    
end


end