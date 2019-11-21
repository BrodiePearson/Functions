function struc = iso_blended_tracer_struc_func_bin(ul,dtheta,r,oceantime,par)
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
if nargin < 5
    par.empty = 'True';
end

dulthetasq = ul.*dtheta.^2;

% Define the bin width
if  isfield(par,'binwid') 
    binwid = par.binwid;
else % default is second order
    binwid = .2;
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
for tt = 1:nt
    for ii = 1:length(bins)

        dulthetasqtmp = dulthetasq(:,tt);
        rtmp = r(:,tt);
        
        if ii == 1
            ind = find(abs(rtmp - bins(ii)) <= abs(rtmp - bins(ii+1)));
        elseif ii == length(bins)
            ind = find(abs(rtmp - bins(ii)) < abs(rtmp - bins(ii-1)));
        else
            ind = find(abs(rtmp - bins(ii)) <= abs(rtmp - bins(ii+1))...
                & abs(rtmp - bins(ii)) < abs(rtmp - bins(ii-1)));   
        end
          
        dulthetasqtmp = dulthetasqtmp(ind);
        
        ulthetasqsftmp(ii) = nanmean(dulthetasqtmp);
        
        bincountstmp(ii) = length(dulthetasqtmp(~isnan(dulthetasqtmp)));
    end

     ulthetasqsf(:,tt) = ulthetasqsftmp;
     
     bincounts(:,tt) = bincountstmp;

    prev=round((tt-1)/nt,2);
    if round(tt/nt,2)~=prev
        disp('...')
    end

end

clear ind
% Take Ensemble Mean along the columns 
% avtsf = nanmean(tsf,2);
% avlsf = nanmean(lsf,2);


sumulthetasq = nansum(ulthetasqsf.*bincounts,2);
totulthetasq = nansum(bincounts,2);
avulthetasq = sumulthetasq./totulthetasq;

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
% [sz,~] = size(tsf);
%  for kk = 1:sz
%     bnct = bincounts(kk,:); 
%     
%     if length(bnct(bnct~=0)) > 1
%         [cil,] = bootci(par.bootsampno,@nanmean,lsf(kk,:));    
%         lci(kk,:) = cil;                       
% 
%         [cit,] = bootci(par.bootsampno,@nanmean,tsf(kk,:));
%         tci(kk,:) = cit ; 
%         
%         ind(kk) = 1;
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
%-------------------------------------------------------------------------
%   Save Data to Structure
%-------------------------------------------------------------------------
% % remove data where there weren't enough timesteps to find a CI
%  ind = logical(ind);
%  avulthetasq = avulthetasq(ind);
%  ulthetasqsf = ulthetasqsf(ind,:);
% %  lci = lci(ind,:);
% %  tci = tci(ind,:);
%  r = bins(ind);
%  
% % remove data where CI is zero for errors with log log plotting
%  ind = tci==0 | lci == 0;
%  ind = sum(ind,2);
%  ind(ind>0) = 1;
%  ind = logical(ind);
%  avulthetasq = avulthetasq(~ind);
%  ulthetasqsf = ulthetasqsf(~ind,:);
% %  lci = lci(~ind,:);
% %  tci = tci(~ind,:);
%  r = r(~ind);
%-------------------------------------------------------------------------
%   Save Data to Structure
%-------------------------------------------------------------------------
struc.bincounts = bincounts;

struc.avulthetasq = avulthetasq;
struc.r = bins;

% struc.lci = lci;
% struc.tci = tci;
 
struc.ulthetasqsf = ulthetasqsf;
struc.oceantime = oceantime;


struc.attributes = ['date modified: ' datestr(date)];

end