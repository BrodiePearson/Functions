function struc = vel_struc_func_bin(ul,ut,r,oceantime,par)
%*************************************************************************
%  struc = vel_struc_func_bin(ul,ut,r,oceantime,par)
%*************************************************************************
%  Function to bin and boostrap pairwise relative velocities and return the
%  isotropic, homogenous, stationary structure function of order n. Output 
%  is a single layer structure.

% np = no. pairs, nt = no. times, nb = no bins
%*************************************************************************
%   Required Input
%-------------------------------------------------------------------------
%	ul,ut       : relative velocities [np,nt] (m/s)
%   r           : pairwise separation distance [np,nt] (km)
%   oceantime   : time [1,nt] (datenumber)
%-------------------------------------------------------------------------
%   Optional Input
%-------------------------------------------------------------------------
%   par         : a structure containing the following variables
%   
%   order       : order of the vel diff, default is 2
%   btstp_time  : boostraps the time mean of the sf
%   btstrp_bins  : boostraps the bins of the structure for each time step
%   bootsampno  : number of bootstrap samples, default is 10,000
%   binwid      : logarithmic bin spacing
%-------------------------------------------------------------------------
%   Default Output
%-------------------------------------------------------------------------
%   lsf,tsf     : longitudinal and transverse structure functions at times
%                  1 through nt [nb,nt] (m/s^order)
%   avlsf/tsf  : time mean of lsf and tsf [nb,1] (m/s^order)
%   r           : separation distance for each pair (in km, if x,y are in
%                 lon/lat)
%   oceantime   : time [1,nt] (datenumber)
%   attributes  : structure containing the date modified
%*************************************************************************
%  Written by Jenna Pearson
%*************************************************************************


%-------------------------------------------------------------------------
%   Setup
%-------------------------------------------------------------------------

[~,nt] = size(ul);

% check to make sure par is defined
if nargin < 6
    par.empty = 'True';
end

% Define the order of the structure function
if  isfield(par,'order') 
    ul = ul.^par.order;
    ut = ut.^par.order;
else % default is second order
    ul = ul.^2;
    ut = ut.^2;
end

% Define the bin width
if  isfield(par,'binwid') 
    binwid = par.binwid;
else % default is second order
    binwid = .2;
    disp('in')
end


%-------------------------------------------------------------------------
%   Bin Vector
%-------------------------------------------------------------------------

nmax = ceil(log10(max(r(:))));
nmin = floor(log10(min(r(:))));
bins = 10.^(nmin:binwid:nmax);

%-------------------------------------------------------------------------
%   Bin Data
%-------------------------------------------------------------------------
disp('Binning data...')

LSF = cell([length(bins),nt]);
TSF = cell([length(bins),nt]);
for tt = 1:nt
    for ii = 1:length(bins)

        ultmp = ul(:,tt);
        uttmp = ut(:,tt);
        rtmp = r(:,tt);
        
        if ii == 1
            ind = find(abs(rtmp - bins(ii)) <= abs(rtmp - bins(ii+1)));
        elseif ii == length(bins)
            ind = find(abs(rtmp - bins(ii)) < abs(rtmp - bins(ii-1)));
        else
            ind = find(abs(rtmp - bins(ii)) <= abs(rtmp - bins(ii+1))...
                & abs(rtmp - bins(ii)) < abs(rtmp - bins(ii-1)));   
        end
        
        ultmp = ultmp(ind);                            
        uttmp = uttmp(ind); 
        
        LSF{ii,tt} = ultmp; 
        TSF{ii,tt} = uttmp; 
        
        lsftmp(ii) = mean(ultmp(~isnan(ultmp)));  
        tsftmp(ii) = mean(uttmp(~isnan(uttmp)));
        
        bincountstmp(ii) = length(uttmp(~isnan(uttmp)));
    end

     tsf(:,tt) = tsftmp;
     lsf(:,tt) = lsftmp;
     bincounts(:,tt) = bincountstmp;

end
 

clear ind

% remove bins where number of pairs are less than 10
ind = bincounts <30;
tsf(ind) = NaN;
lsf(ind) = NaN;
r(ind) = NaN;
bincounts(ind) = NaN;

clear ind

% Take Ensemble Mean along the columns 
avtsf = nanmean(tsf,2);
avlsf = nanmean(lsf,2);

% sumtsf = nansum(tsf.*bincounts,2);
% tottsf = nansum(bincounts,2);
% avtsf = sumtsf./tottsf;
% 
% sumlsf = nansum(lsf.*bincounts,2);
% totlsf = nansum(bincounts,2);
% avlsf = sumlsf./totlsf;

%-------------------------------------------------------------------------
%   Bootstrapping
%-------------------------------------------------------------------------
disp('Boostrapping...')

% default number of resamples is 10000
if  ~isfield(par,'bootsampno')
    par.bootsampno = 10000;
end

[sz,nt] = size(tsf);

if nt >1
    lci = NaN([sz,2]);
    tci = NaN([sz,2]);
    for kk = 1:sz
        bnct = bincounts(kk,:); 

        if nansum(bnct) > 30 % This corresponds to have at least 8 drifters or locations
            [cil,] = bootci(par.bootsampno,@nanmean,lsf(kk,:));    
            lci(kk,:) = cil;                       

            [cit,] = bootci(par.bootsampno,@nanmean,tsf(kk,:));
            tci(kk,:) = cit ; 
        end
     end
else
    for tt = 1:nt
        
        for bb =  1:length(bins)
            bnct = bincounts(bb,tt); 
            if bnct > 30 && ~isnan(bnct) % This corresponds to have at least 8 drifters or locations
                [cil,] = bootci(par.bootsampno,@nanmean,LSF{bb,tt});    
                LCI(bb,:) = cil;                       

                [cit,] = bootci(par.bootsampno,@nanmean,TSF{bb,tt});
                TCI(bb,:) = cit ; 
            else
                LCI(bb,:) = NaN([1,2]);
                TCI(bb,:) = NaN([1,2]);
            end
            
        end
     end
end
%-------------------------------------------------------------------------
%   Save Data to Structure
%-------------------------------------------------------------------------

if nt >1
    % remove data where CI is NaN
    ind = isnan(tci) | isnan(lci);
    ind = sum(ind,2);
    ind(ind>0) = 1;
    ind = logical(ind);
    avtsf = avtsf(~ind);
    avlsf = avlsf(~ind);
    tsf = tsf(~ind,:);
    lsf = lsf(~ind,:);
    lci = lci(~ind,:);
    tci = tci(~ind,:);
    r = bins(~ind);
    bincounts = bincounts(~ind,:);
else
    % remove data where CI is NaN
    ind = isnan(TCI) | isnan(LCI);
    ind = sum(ind,2);
    ind(ind>0) = 1;
    ind = logical(ind);
    avtsf = avtsf(~ind);
    avlsf = avlsf(~ind);
    tsf = tsf(~ind,:);
    lsf = lsf(~ind,:);
    LCI = LCI(~ind,:);
    TCI = TCI(~ind,:);
    r = bins(~ind);
    bincounts = bincounts(~ind,:);
end


%-------------------------------------------------------------------------
%   Save Data to Structure
%-------------------------------------------------------------------------
disp('Saving Data...')
struc.bincounts = bincounts;

struc.avtsf = avtsf;
struc.avlsf = avlsf;

struc.tsf = tsf;
struc.lsf = lsf; 

struc.LSF = LSF;
struc.TSF = TSF;

struc.r = r;
% struc.R = R;

struc.oceantime = oceantime;

if nt >1
    struc.lci = lci;
    struc.tci = tci;
else
    struc.LCI = LCI;
    struc.TCI = TCI;
end

struc.attributes = {['date modified: ' datestr(date)];'Author: Jenna Pearson'; 'Contact: jenna_pearson@brown.edu'};
struc.par = par;

end