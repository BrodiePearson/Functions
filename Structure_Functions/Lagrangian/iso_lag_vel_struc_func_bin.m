function struc = iso_lag_vel_struc_func_bin(du,dv,dul,dut,tau,par)
%*************************************************************************
%  struc = lvel_struc_func_bin(ul,ut,r,oceantime,par)
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
%   bootsampno  : number of bootstrap samples, default is 10,000
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

[~,nl] = size(du);

% check to make sure par is defined
if nargin < 6
    par.empty = 'True';
end

% Define the order of the structure function
if  isfield(par,'order') 
    du = du.^par.order;
    dv = dv.^par.order;
    dul = dul.^par.order;
    dut = dut.^par.order;
else % default is second order
    du = du.^2;
    dv = dv.^2;
    dul = dul.^2;
    dut = dut.^2;
end

if isfield(par,'binwid') 
    binwid = par.binwid;
else
    binwid = 0.1; % default
end

%-------------------------------------------------------------------------
%   Bin Vector
%-------------------------------------------------------------------------

nmax = ceil(log10(max(tau(:))));
nmin = floor(log10(min(tau(:))));
bins = 10.^(nmin:binwid:nmax);

%-------------------------------------------------------------------------
%   Bin Data
%-------------------------------------------------------------------------
body('Binning data.')
for ll = 1:nl
    for ii = 1:length(bins)

        dutmp = du(:,ll);
        dvtmp = dv(:,ll);
        dultmp = dul(:,ll);
        duttmp = dut(:,ll);
        tautmp = tau(:,ll);
        
        if ii == 1
            ind = find(abs(tautmp - bins(ii)) <= abs(tautmp - bins(ii+1)));
        elseif ii == length(bins)
            ind = find(abs(tautmp - bins(ii)) < abs(tautmp - bins(ii-1)));
        else
            ind = find(abs(tautmp - bins(ii)) <= abs(tautmp - bins(ii+1))...
                & abs(tautmp - bins(ii)) < abs(tautmp - bins(ii-1)));   
        end
        
        dutmp = dutmp(ind);                            
        dvtmp = dvtmp(ind);  
        dultmp = dultmp(ind);
        duttmp = duttmp(ind);
        
        lusftmp(ii) = mean(dutmp(~isnan(dutmp)));  
        lvsftmp(ii) = mean(dvtmp(~isnan(dvtmp)));
        lulsftmp(ii) = mean(dultmp(~isnan(dultmp)));  
        lutsftmp(ii) = mean(duttmp(~isnan(duttmp)));
        
        bincountstmp(ii) = length(dvtmp(~isnan(dvtmp)));
    end

    usf(:,ll) = lusftmp;
    vsf(:,ll) = lvsftmp;
    ulsf(:,ll) = lulsftmp;
    utsf(:,ll) = lutsftmp;
    bincounts(:,ll) = bincountstmp;

end

% Take Ensemble Mean along the columns 
avusf = nanmean(usf,2);
avvsf = nanmean(vsf,2);
avulsf = nanmean(ulsf,2);
avutsf = nanmean(utsf,2);
%-------------------------------------------------------------------------
%   Bootstrapping
%-------------------------------------------------------------------------
body('Boostrapping.')

% default number of resamples is 10000
if  ~isfield(par,'bootsampno')
    par.bootsampno = 10000;
end

[ntdiff,ntraj] = size(usf);

if ntraj >1
    uci = NaN([ntdiff,2]);
    vci = NaN([ntdiff,2]);
    ulci = NaN([ntdiff,2]);
    utci = NaN([ntdiff,2]);
    
    for kk = 1:ntdiff
        bnct = bincounts(kk,:); 

        if nansum(bnct) > 30 % This corresponds to at least 8 time differences
            [ciu,] = bootci(par.bootsampno,@nanmean,usf(kk,:));    
            uci(kk,:) = ciu;                       

            [civ,] = bootci(par.bootsampno,@nanmean,vsf(kk,:));
            vci(kk,:) = civ ; 
            
            [ciul,] = bootci(par.bootsampno,@nanmean,ulsf(kk,:));    
            ulci(kk,:) = ciul;   
            
            [ciut,] = bootci(par.bootsampno,@nanmean,utsf(kk,:));    
            utci(kk,:) = ciut;   
        end
     end
else
    for dd = 1:ntraj
        
        for bb =  1:length(bins)
            bnct = bincounts(bb,dd); 
            if bnct > 30 && ~isnan(bnct) % This corresponds to have at least 8 timesteps
                [ciu,] = bootci(par.bootsampno,@nanmean,USF{bb,dd});    
                UCI(bb,:) = ciu;                       

                [civ,] = bootci(par.bootsampno,@nanmean,VSF{bb,dd});
                VCI(bb,:) = civ ; 
                
                [ciul,] = bootci(par.bootsampno,@nanmean,ULSF{bb,dd});    
                ULCI(bb,:) = ciul; 
                
                [ciut,] = bootci(par.bootsampno,@nanmean,UTSF{bb,dd});    
                UTCI(bb,:) = ciut; 
            else
                UCI(bb,:) = NaN([1,2]);
                VCI(bb,:) = NaN([1,2]);
                ULCI(bb,:) = NaN([1,2]);
                UTCI(bb,:) = NaN([1,2]);
            end
            
        end
     end
end
%-------------------------------------------------------------------------
%   Get rid of bad data
%-------------------------------------------------------------------------

if ntraj >1
    % remove data where CI is NaN
    ind = isnan(vci) | isnan(uci) | isnan(ulci) | isnan(utci);
    ind = sum(ind,2);
    ind(ind>0) = 1;
    ind = logical(ind);
    avusf = avusf(~ind);
    avvsf = avvsf(~ind);
    avulsf = avulsf(~ind);
    avutsf = avutsf(~ind);
    usf = usf(~ind,:);
    vsf = vsf(~ind,:);
    ulsf = ulsf(~ind,:);
    utsf = utsf(~ind,:);
    uci = uci(~ind,:);
    vci = vci(~ind,:);
    ulci = ulci(~ind,:);
    utci = utci(~ind,:);
    tau = bins(~ind);
    bincounts = bincounts(~ind,:);
else
    % remove data where CI is NaN
    ind = isnan(VCI) | isnan(UCI);
    ind = sum(ind,2);
    ind(ind>0) = 1;
    ind = logical(ind);
    avusf = avusf(~ind);
    avvsf = avvsf(~ind);
    avulsf = avulsf(~ind);
    avutsf = avutsf(~ind);
    usf = usf(~ind,:);
    vsf = vsf(~ind,:);
    ulsf = ulsf(~ind,:);
    ulsf = ulsf(~ind,:);
    UCI = UCI(~ind,:);
    VCI = VCI(~ind,:);
    ULCI = ULCI(~ind,:);
    UTCI = UTCI(~ind,:);
    tau = bins(~ind);
    bincounts = bincounts(~ind,:);
end

%-------------------------------------------------------------------------
%   Save Data to Structure
%-------------------------------------------------------------------------
body('Saving data to structure.')
struc.par = par;
struc.bincounts = bincounts;

struc.avusf = avusf;
struc.avvsf = avvsf;
struc.avulsf = avulsf;
struc.avutsf = avutsf;
struc.tau = tau;

struc.uci = uci;
struc.vci = vci;
struc.ulci = ulci;
struc.utci = utci;

struc.usf = usf;
struc.vsf = vsf;
struc.ulsf = ulsf; 
struc.utsf = utsf;

struc.attributes = ['date modified: ' datestr(date)];


% %-------------------------------------------------------------------------
% %   Bootstrapping
% %-------------------------------------------------------------------------
% disp('Boostrapping...')
% 
% % default number of resamples is 10000
% if  ~isfield(par,'bootsampno')
%     par.bootsampno = 10000;
% end
% 
% [sz,~] = size(vsf);
%  for kk = 1:sz
%     bnct = bincounts(kk,:); 
%     
%     if length(bnct(bnct~=0)) > 1
%         [cil,] = bootci(par.bootsampno,@nanmean,usf(kk,:));    
%         uci(kk,:) = cil;                       
% 
%         [cit,] = bootci(par.bootsampno,@nanmean,vsf(kk,:));
%         vci(kk,:) = cit ; 
%         
%         ind(kk) = 1;
%         
%     else
%         ind(kk) = 0;
%     end
%     
%     prev=round((kk-1)/sz,2);
% 
%     if round(kk/sz,2)~=prev
%         disp([int2str(round(kk/sz,2)*100) '%'])
%     end
%  end
%  
% %-------------------------------------------------------------------------
% %   Save Data to Structure
% %-------------------------------------------------------------------------
% % remove data where there weren't enough timesteps to find a CI
%  ind = logical(ind);
%  avvsf = avvsf(ind);
%  avusf = avusf(ind);
%  uci = uci(ind,:);
%  vci = vci(ind,:);
%  tau = bins(ind);
%  
% % remove data where CI is zero for errors with log log plotting
%  ind = vci==0 | uci == 0;
%  ind = sum(ind,2);
%  ind(ind>0) = 1;
%  ind = logical(ind);
%  avvsf = avvsf(~ind);
%  avusf = avusf(~ind);
%  uci = uci(~ind,:);
%  vci = vci(~ind,:);
%  tau= tau(~ind);
 
end