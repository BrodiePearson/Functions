function struc = tracer_struc_func_bin(dbeta,dx,dy,oceantime,par)
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

[~,nt] = size(dbeta);

% check to make sure par is defined
if nargin < 5
    par.empty = 'True';
end

dbetasq = dbeta.^2;

% Define the bin width
if  isfield(par,'binwid') 
    if par.binwid < min(dx(:)) || par.binwid < min(dy(:))
        error('Error: binwidth to small.')
    else
        binwid = par.binwid;
    end
else % default 1
    binwid = 1;
    par.binwid = binwid;
end


%-------------------------------------------------------------------------
%   Bin Vector
%-------------------------------------------------------------------------
body('Generating bin vectors.')

nmax = ceil(max(dx(:)));
nmin = floor(min(dx(:)));
binsx = nmin:binwid:nmax;

nmax = ceil(max(dy(:)));
nmin = floor(min(dy(:)));
binsy = nmin:binwid:nmax;

%-------------------------------------------------------------------------
%   Initialize matrices
%-------------------------------------------------------------------------
betasqsf = NaN([length(binsx),length(binsy),size(dbetasq,3)]);
bincounts = NaN([length(binsx),length(binsy),size(dbetasq,3)]);
%-------------------------------------------------------------------------
%   Bin Data
%-------------------------------------------------------------------------
body('Binning data.')
for tt = 1:nt
    for xx = 1:length(binsx)
        for yy = 1:length(binsy)
            dbetasqtmp = dbetasq(:,tt);
            dxtmp = dx(:,tt);
            dytmp = dy(:,tt); 

            % upper left corner
            if xx == 1 && yy == 1
                ind = find(abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx+1)) & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy+1)));
                
            % lower right corner
            elseif xx == length(binsx) && yy == length(binsy)
                ind = find(abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx-1)) & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy-1)));
                
            % lower left corner
            elseif xx == 1 && yy == length(binsy)
                ind = find(abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx+1))...
                   & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy-1)));
            % upper right corner
            elseif yy == 1 && xx == length(binsx)
                ind = find(abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx-1))...
                    & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy+1)));
            % left colmun   
            elseif xx == 1
                ind = find(abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx+1))...
                    & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy+1))...
                    & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy-1))); 
            % top row
            elseif yy == 1
                ind = find(abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx+1))...
                    & abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx-1))...
                    & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy+1)));
                
            % right column
            elseif xx == length(binsx)
                ind = find(abs(dxtmp - binsx(xx)) < abs(dxtmp - binsx(xx-1))...
                    & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy+1))...
                    & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy-1)));
            % bottom row
            elseif yy == length(binsy)
                ind = find(abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx+1))...
                    & abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx-1))...
                    & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy-1)));
            % interior points
            else
                ind = find(abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx+1))...
                    & abs(dxtmp - binsx(xx)) <= abs(dxtmp - binsx(xx-1))...
                    & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy+1))...
                    & abs(dytmp - binsy(yy)) <= abs(dytmp - binsy(yy-1)));   
            end

            dbetasqtmp = dbetasqtmp(ind);
            betasqsftmp(xx,yy) = nanmean(dbetasqtmp);
            bincountstmp(xx,yy) = length(dbetasqtmp(~isnan(dbetasqtmp)));
        end
    end
    
    betasqsf(:,:,tt) = betasqsftmp;
    bincounts(:,:,tt) = bincountstmp;

end

% Take Ensemble Mean along the columns 
avbetasqsf = nanmean(betasqsf,3);

%-------------------------------------------------------------------------
%   Save Data to Structure
%-------------------------------------------------------------------------
body('Saving data to structure.')
struc.bincounts = bincounts;

struc.avbetasqsf = avbetasqsf';
struc.binsx = binsx;
struc.binsy = binsy;
struc.betasqsf = betasqsf;
struc.oceantime = oceantime;
struc.par = par;
end