function struc = blended_tracer_struc_func_bin(dx,dy,ul,ut,dbeta,oceantime,par)
%*************************************************************************
%  Function to bin and boostrap pairwise relative velocities and return the
%  isotropic, homogenous, stationary structure function of order n. Output 
%  is a single layer structure.

% np = no. pairs, nt = no. times, nb = no bins
%*************************************************************************
%  Syntax:
%-------------------------------------------------------------------------
%   struc = blended_tracer_struc_func_bin(dx,dy,ul,ut,dbeta,oceantime,Property,Value,...)
%*************************************************************************
%   Required Input
%-------------------------------------------------------------------------
%	ul,ut       : relative velocities [np,nt] (m/s)
%   r           : pairwise separation distance [np,nt] (km)
%   oceantime   : time [1,nt] (datenumber)
%-------------------------------------------------------------------------
%   Optional Input
%-------------------------------------------------------------------------   
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
% Error Checks and Parse Optional Inputs
%-------------------------------------------------------------------------
% narginchk(6,inf)
% assert(isequal(size(x),size(y),size(u),size(v),size(beta))==1,'Error: Dimensions of x, y, u, v, and beta must agree. Check inputs, or use [x,y] = meshgrid(x,y) if necessary.')
% assert(isequal(size(x,2),length(oceantime)),'Error: Length of variable oceantime must match second dimension of x, y, u, v, and beta.') 
% 
% tmp = strncmpi(varargin,'inhomogeneity',3); 
% if any(tmp)
%    inhomogeneity = varargin{find(tmp)+1}; 
%    tmp(find(tmp)+1)=1; 
%    varargin = varargin(~tmp); 
%    assert(islogical(inhomogeneity),'Input error: Inhomogeneity option must be true or false. Default is true.') 
% end
% 
% tmp = strncmpi(varargin,'anisotropy',3); 
% if any(tmp)
%    anisotropy = varargin{find(tmp)+1}; 
%    tmp(find(tmp)+1)=1; 
%    varargin = varargin(~tmp); 
%    assert(islogical(anisotropy),'Input error: Anisotropy option must be true or false. Default is true.') 
% end

%-------------------------------------------------------------------------
%   Setup
%-------------------------------------------------------------------------

[~,nt] = size(ul);
 
dulbetasq = ul.*dbeta.^2;
dutbetasq = ut.*dbeta.^2;

% Define the bin width
if  isfield(par,'binwid') 
    binwid = par.binwid;
else % default is second order
    binwid = 1;
end


%-------------------------------------------------------------------------
%   Bin Vector
%-------------------------------------------------------------------------

% nmax = ceil(log10(abs(max(dx(:)))));
% nmin = floor(log10(abs(min(dx(:)))));
% binsx = 10.^(-nmin:binwid:nmax);
% 
% nmax = ceil(log10(abs(max(dy(:)))));
% nmin = floor(log10(abs(min(dy(:)))));
% binsy = 10.^(-nmin:binwid:nmax);

nmax = ceil(max(dx(:)));
nmin = floor(min(dx(:)));
binsx = nmin:binwid:nmax;

nmax = ceil(max(dy(:)));
nmin = floor(min(dy(:)));
binsy = nmin:binwid:nmax;


%-------------------------------------------------------------------------
%   Initialize matrices
%-------------------------------------------------------------------------
ulbetasqsf = NaN([length(binsx),length(binsy),size(dulbetasq,3)]);
utbetasqsf = NaN([length(binsx),length(binsy),size(dutbetasq,3)]);
bincounts = NaN([length(binsx),length(binsy),size(dulbetasq,3)]);


%-------------------------------------------------------------------------
%   Bin Data
%-------------------------------------------------------------------------
disp('Binning data...')
parfor tt = 1:nt
    
    ulbetasqsftmp = NaN([length(binsx),length(binsy)]);
    utbetasqsftmp = NaN([length(binsx),length(binsy)]);
    bincountstmp = NaN([length(binsx),length(binsy)]);

    for xx = 1:length(binsx)
        for yy = 1:length(binsy)
            dulbetasqtmp = dulbetasq(:,tt);
            dutbetasqtmp = dutbetasq(:,tt);
            dxtmp = dx(:,tt);
            dytmp = dy(:,tt); 

            
            if xx == 1 && yy == 1
                ind = find(abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx+1)) & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy+1)));
                
                
            elseif xx == length(binsx) && yy == length(binsy)
                ind = find(abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx-1)) & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy-1)));
                
            elseif xx == 1 && yy == length(binsy)
                ind = find(abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx+1))...
                   & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy-1)));
               
              
            elseif yy == 1 && xx == length(binsx)
                ind = find(abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx-1))...
                    & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy+1)));
                
            elseif xx == 1
                ind = find(abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx+1))...
                    & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy+1))...
                    & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy-1))); 
            elseif yy == 1
                ind = find(abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx+1))...
                    & abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx-1))...
                    & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy+1)));
                
                
            elseif xx == length(binsx)
                ind = find(abs(dxtmp - binsx(xx)) < abs(dxtmp - binsx(xx-1))...
                    & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy+1))...
                    & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy-1)));
            elseif yy == length(binsy)
                ind = find(abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx+1))...
                    & abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx-1))...
                    & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy-1)));
                
                
            else
                ind = find(abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx+1))...
                    & abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx-1))...
                    & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy+1))...
                    & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy-1)));   
            end

            dulbetasqtmp = dulbetasqtmp(ind);
            dutbetasqtmp = dutbetasqtmp(ind);
            
            ulbetasqsftmp(xx,yy) = nanmean(dulbetasqtmp);
            utbetasqsftmp(xx,yy) = nanmean(dutbetasqtmp);
            
            bincountstmp(xx,yy) = length(dulbetasqtmp(~isnan(dulbetasqtmp)));
        end
    end
    
    ulbetasqsf(:,:,tt) = ulbetasqsftmp;
    utbetasqsf(:,:,tt) = utbetasqsftmp;
    bincounts(:,:,tt) = bincountstmp;

end

clear ind

% Take Ensemble Mean along the columns 
avulbetasqsf = nanmean(ulbetasqsf,3);
avutbetasqsf = nanmean(utbetasqsf,3);

% sumulbetasq = nansum(ulbetasqsf.*bincounts,3);
% totulbetasq = nansum(bincounts,3);
% avulbetasq = sumulbetasq./totulbetasq;

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

struc.avulbetasqsf = avulbetasqsf';
struc.avutbetasqsf = avutbetasqsf';
struc.binsx = binsx;
struc.binsy = binsy;

struc.binwid = binwid;
 
struc.ulbetasqsf = ulbetasqsf;
struc.utbetasqsf = utbetasqsf;
struc.oceantime = oceantime;

end