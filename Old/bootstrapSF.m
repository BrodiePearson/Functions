function SFboot = bootstrapSF(r, TSF,LSF,par)
%% Info

% USE:
%       To bootstrap outputs of means from VSF_Calc
% INPUT:
%   TSF = nd x nt longitudinal structure function values (m/s)^2
%   LSF = nd x nt transverse structure function values (m/s)^2
%   par = structure of parameters
%       par.n = integer number of bootstrap samples (e.g. 10000)

% OUTPUT:

%% Bootstrapping
disp('Boostrapping...')

[sz,nt] = size(TSF);

if par.interpolate == 1
     for ii = 1:nt
        [lsfc,gof,output] = fit(r,LSF(:,ii),'smoothingspline');
        [tsfc,gof,output] = fit(r,TSF(:,ii),'smoothingspline');


        x = (min(r):max(r))';
        yl = feval(lsfc,x);
        yt = feval(tsfc,x);

        yl(1) = LSF(1,ii);
        yt(1) = TSF(1,ii);

        yl(end) = LSF(end,ii);
        yt(end) = TSF(end,ii);
        tempLSF(:,ii) = yl;
        tempTSF(:,ii) = yt;
        tempr(:,ii) = x;
     end
     
     LSF = tempLSF;
     TSF = tempTSF;
     r = tempr;
     [sz,nt] = size(TSF);
end

for kk = 1:sz
 
    [btcilong,] = bootci(par.n,@nanmean,LSF(kk,:));    % Bootstrap LSF to find confidence intervals
    lci(kk,1) = btcilong(1,1);                       % Save lower confidence interval value
    lci(kk,2) = btcilong(2,1);                       % Save upper confidence interval value

    [btcitran,] = bootci(par.n,@nanmean,TSF(kk,:));    % Bootstrap TSF to find confidence intervals
    tci(kk,1) = btcitran(1,1);                       % Save lower confidence interval value
    tci(kk,2) = btcitran(2,1);                       % Save upper confidence interval valu
    prev=round((kk-1)/sz,2);
    
    if round(kk/sz,2)~=prev
        disp(['PROGRESS: ', int2str(round(kk/sz,2)*100), '%'])
    end
end

SFboot.l = LSF;
SFboot.t = TSF;
SFboot.r = r;
SFboot.lCIL= lci(:,1);
SFboot.lCIU= lci(:,2);
SFboot.tCIL = tci(:,1);
SFboot.tCIU = tci(:,2); 


end