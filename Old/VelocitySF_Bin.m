function VSF_binned = VelocitySF_Bin(l,t,r,oceantime,par)
%% Info:

% USE:

%   To bin and bootstrap the calculated isotropic, homogenous, stationary 
%   velocity structure function output from VelocitySF_Calc

% INPUT: nd x nt = matrix of dimensions number velocity differences x number of time steps
%       l = nd x nt longitudinal velocity differences (m/s)
%       t = nd x nt transverse velocity differences (m/s)
%       r = nd x nt separation distances (km)
%       oceantime = 1 x nt row vector of time steps (Datenumber)
%       par = structure with set of parameters 
%           par.binexp = intial exponent of bins [larger = smaller bins represented, recommended 10] (integer)
%           par.order = order of structure function
%           par.bootsampnum = # of bootstrap resamples to choose [recommended 10,000](integer)

% OPTIONAL INPUT:
%           par.subsample = 'True' will return colindx and rowindx from VelocitySF_calc ('True'/'False)

% OUTPUT: nb = number of bins
%       l = 1 x nb longitudinal structure function values (m/s)^n
%       t = 1 x nb transverse structure function values (m/s)^n
%       r = 1 x nb logarithmically spaced, centered separation distances (km)
%       oceantime = 1 x nt row vector of time steps (Datenumber)
%       TSF = nb x nt matrix of tsf values in each bin before averaging
%       LSF = nb x nt matrix of lsf values in each bin before averaging
%       lCIL/tCIL = 1 x nb lower bootstrapped confidence interval estimate 
%       lCIU/tCIU = 1 x nb upper bootstrapped confidence interval estimate
%       bincounts = nb x nt matrix with the number of points in each bin at every time step

% OPTIONAL OUTPUT:
%       If par.subsample = 'True'
%           colindx,rowindx = subsampling indices from VelocitySF_calc  

%% Setup

% Define the order of the structure function
l = l.^par.order;
t = t.^par.order;

%% Structure Function Binning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create a logarithmically spaced bin vector that increases by a power of 10 every four points
%% Structure Function Binning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Bins(1) = 0;
kk=2;
while(10^((kk-par.binexp)*.25) <= ceil(max(r(:))))
    Bins(kk) = 10^((kk-par.binexp)*.25); % Logarithmically spaced bin vector that increases by a power of 10 every four points
    kk = kk+1;
end

% Time-stamped structure functions
for tt = 1:length(oceantime)

    for ii = 1:length(Bins)-1

        ind = find(r(:,tt) < Bins(ii+1) & r(:,tt) >= Bins(ii));   % Isolate the indices of data that fall into a given bin

        longvel = l(:,tt);
        tranvel = t(:,tt);

        sflv = longvel(ind);                            % Find long. sf values associated with the indices from ^
        sftv = tranvel(ind);                            % Find tran. sf values associated with the indices from ^
        tsf(ii) = mean(sftv(~isnan(sftv)));                        % Compute the average of all data in the tran. sf bin and save to vector for plotting
        lsf(ii) = mean(sflv(~isnan(sflv)));                        % Compute the average of all data in the long. sf bin and save to vector for plotting
        counttsf(ii) = length(sftv(~isnan(sftv)));
        countlsf(ii) = length(sflv(~isnan(sflv)));

    end

     % Add info to matrix

     TSF(:,tt) = tsf;
     LSF(:,tt) = lsf;
     COUNTTSF(:,tt) = counttsf;
     COUNTLSF(:,tt) = countlsf;

        prev=round((tt-1)/length(oceantime),2);
        if round(tt/length(oceantime),2)~=prev
            disp(['PROGRESS: ', int2str(round(tt/length(oceantime),2)*100), '%'])
        end

end

% Take Ensemble Mean along the columns 

% sumtsf = nansum(TSF.*COUNTTSF,2);
% tottsf = nansum(COUNTTSF,2);
% enstsf = sumtsf./tottsf;
% 
% sumlsf = nansum(LSF.*COUNTLSF,2);
% totlsf = nansum(COUNTLSF,2);
% enslsf = sumlsf./totlsf;

enstsf = nanmean(TSF,2);
enslsf = nanmean(LSF,2);

% Center Bins
disp('Centering Bins...')
for jj = 1:length(Bins)-1
    rmid(jj) = (Bins(jj)+Bins(jj+1))/2;     % Declare the average distance between bins as the representative displacement for a given bin ( centering )
end

%% Bootstrapping

disp('Boostrapping...')

[sz,~] = size(TSF);

 for kk = 1:sz
    if isnan(nanmean(LSF(kk,:))) == 0
        [btcilong,] = bootci(par.bootsampnum,@nanmean,LSF(kk,:));    % Bootstrap LSF to find confidence intervals
        lci(kk,1) = btcilong(1,1);                       % Save lower confidence interval value
        lci(kk,2) = btcilong(2,1);                       % Save upper confidence interval value

        [btcitran,] = bootci(par.bootsampnum,@nanmean,TSF(kk,:));    % Bootstrap TSF to find confidence intervals
        tci(kk,1) = btcitran(1,1);                       % Save lower confidence interval value
        tci(kk,2) = btcitran(2,1);                       % Save upper confidence interval valu
    else
        lci(kk,1) = NaN;                       % Save lower confidence interval value
        lci(kk,2) = NaN;                       % Save upper confidence interval value
        tci(kk,1) = NaN;                       % Save lower confidence interval value
        tci(kk,2) = NaN;                       % Save upper confidence interval valu
    end
        prev=round((kk-1)/sz,2);

    if round(kk/sz,2)~=prev
        disp([int2str(round(kk/sz,2)*100) '%'])
    end
end





%% Write Data
VSF_binned.bincounts = COUNTTSF;

VSF_binned.TSF = TSF;
VSF_binned.LSF = LSF;
VSF_binned.R = r;  

VSF_binned.oceantime = oceantime;

VSF_binned.lCIL= lci(:,1);
VSF_binned.lCIL = VSF_binned.lCIL(~isnan(enslsf));
VSF_binned.lCIU = lci(:,2);
VSF_binned.lCIU = VSF_binned.lCIU(~isnan(enslsf));
VSF_binned.tCIL = tci(:,1);
VSF_binned.tCIL = VSF_binned.tCIL(~isnan(enslsf));
VSF_binned.tCIU = tci(:,2); 
VSF_binned.tCIU = VSF_binned.tCIU(~isnan(enslsf));

VSF_binned.l = enslsf(~isnan(enslsf));
VSF_binned.t = enstsf(~isnan(enslsf));
VSF_binned.r = rmid(~isnan(enslsf))';
VSF_binned.attributes = ['date modified: ' datestr(date)];

if strcmp(par.subsample,'True') == 1
    VSF_binned.colindx=par.colindx;
    VSF_binned.rowindx=par.rowindx;
end

end