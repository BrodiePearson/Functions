function VSF2D_binnned = VSF2D_Bin(l,t,ds1,ds2,oceantime,par)
%% Info:

% USE:

%   To calculate the angle-dependent velocity structure function of any order that is of the form
%   [delta(u)^n]

% INPUT: nd x nt = matrix of dimensions number velocity differences x number of time steps
%       ds1 = nd X nt, spatial difference matrix (dx or dr)
%       ds2 = nd x nt, spatial difference matrix (dy or dtheta)
%       dl = nd x nt longitudinal velocity difference matix       
%       dt = nd x nt transverse velocity difference matix 
%       oceantime = row timevector (1 x nt)
%       par = structure with set of parameters 
%           par.order = order of structure function
%           par.spacesample = percent of nd to subsample
%           par.timesample = percent of nt to subsample
%           par.resamp = percent of subsampled data to pad
%           par.padding = number of grid points to pad for each resamp pt
%           par.direction = 
%               'xy' -> zonal.meridional
%               'rt' -> r, theta (coming soon)
%           par.binint = interval between linear bins (in km)

% OUTPUT: 


%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

order = par.order;
spacesample = par.spacesample;
timesample = par.timesample;
percentresamp = par.percentresamp;
padding = par.padding;
direction = par.direction;
binint = par.binint;


% subsample time
[ ~,colindx] = datasample(oceantime,timesample,2,'Replace',false); % datasample(dataset,number of samples, dimension [1=row,2=col],...)
colindx =sort(colindx);

ds1 = ds1(:,colindx);               
ds2 = ds2(:,colindx);
dv1 = dv1(:,colindx);
dv2 = dv2(:,colindx);
oceantime = oceantime(colindx);

clear colindx

% Randomly subsample

[ ~,rowindx] = datasample(ds1,spacesample,1,'Replace',false); % datasample(dataset,number of samples, dimension [1=row,2=col],...)


resamp = round(percentresamp/100*length(rowindx),2);
ind = datasample(rowindx,resamp,'Replace',false);

% add in grid points plus or minus 4 grid points away from those
% originally sampled

for aa=1:padding
    indp = ind + aa.*ones(1,resamp);
    indn = ind - aa.*ones(1,resamp);
    rowindx = [rowindx indp indn];
end


clear ind indp indn resamp

rowindx =sort(rowindx);

ds1 = ds1(rowindx,:);               
ds2 = ds2(rowindx,:);
dv1 = dv1(rowindx,:);
dv2 = dv2(rowindx,:);

clear rowindx

[nd,nt] = size(ds1);


%% Calculations
if strcmp(direction,'xy') == 1

    dx = ds1;
    dy = ds2;

    Binsx = floor(min(dx(:))):binint:ceil(max(dx(:)));   
    Binsy = floor(min(dy(:))):binint:ceil(max(dy(:)));   

    for ii = 1:length(Binsx)-1
        for jj = 1:length(Binsy)-1

            ind = find(dx <= Binsx(ii+1) & dx > Binsx(ii) & dy <= Binsy(jj+1) & dy > Binsy(jj));   % Isolate the indices of data that fall into a given bin
            sflvxy = l(ind);                            % Find long. sf values associated with the indices from ^
            sftvxy = t(ind);                            % Find tran. sf values associated with the indices from ^
            tsfxy(ii,jj) = mean(sftvxy);                      % Compute the average of all data in the tran. sf bin and save to vector for plotting
            lsfxy(ii,jj) = mean(sflvxy);                      % Compute the average of all data in the long. sf bin and save to vector for plotting
            counttsfxy(ii,jj) = length(sftvxy);               % Count number of elements in each bin 
            countlsfxy(ii,jj) = length(sflvxy);
            clear sflvxy sftvxy ind
        end

        prev=round((ii-1)/(length(Binsx)-1),2);
        if round(ii/(length(Binsx)-1),2)~=prev
            disp(['PROGRESS: ', int2str(round(ii/(length(Binsx)-1),2)*100), '%'])
        end

    end


    disp('Centering Bins...')
    for jj = 1:length(Binsx)-1
        rmidx(jj) = (Binsx(jj)+Binsx(jj+1))/2;     % Declare the average distance between bins as the representative displacement for a given bin ( centering )
    end

    for jj = 1:length(Binsy)-1
        rmidy(jj) = (Binsy(jj)+Binsy(jj+1))/2;     % Declare the average distance between bins as the representative displacement for a given bin ( centering )
    end

    % Save data to structure
    VSF2D_binned.tsf = tsfxy;
    VSF2D_binned.lsf = lsfxy;
    VSF2D_binned.rdx = rmidx;  
    VSF2D_binned.rdy = rmidy;

    VSF2D_binned.counttsf = counttsfxy;
    VSF2D_binned.countlsf = countlsfxy;

else
    Binsr = min(min(dr)):.1:max(max(dr));             % 100m bins
    Binstheta = min(min(dtheta)):.0001:max(max(dtheta));   % .001 degree bins 
    
    
    for ii = 1:length(Binsr)-1
        for jj = 1:length(Binstheta)-1
        
            ind = find(dr <= Binsr(ii+1) & dr > Binsr(ii) & dtheta <= Binstheta(jj+1) & dtheta > Binstheta(jj));   % Isolate the indices of data that fall into a given bin
            sflvrtheta = l(ind);                            % Find long. sf values associated with the indices from ^
            sftvrtheta = t(ind);                            % Find tran. sf values associated with the indices from ^
            tsfrtheta(ii,jj) = mean(sftvrtheta);                      % Compute the average of all data in the tran. sf bin and save to vector for plotting
            lsfrtheta(ii,jj) = mean(sflvrtheta);                      % Compute the average of all data in the long. sf bin and save to vector for plotting
            counttsfrtheta(ii,jj) = length(sftvrtheta);               % Count number of elements in each bin 
            countlsfrtheta(ii,jj) = length(sflvrtheta);
        end
        
        prev=round((ii-1)/(length(Binsr)-1),2);
        if round(ii/(length(Binsr)-1),2)~=prev
            disp(['PROGRESS FOR BINNING SF VALUES: ', int2str(round(ii/(length(Binsr)-1),2)*100), '%'])
        end
        
        clear sflvrtheta sftvrtheta
    end
    
    
    disp('Centering Bins...')
    for jj = 1:length(Binsr)-1
        rmidr(jj) = (Binsr(jj)+Binsr(jj+1))/2;     % Declare the average distance between bins as the representative displacement for a given bin ( centering )
    end
    
    for jj = 1:length(Binstheta)-1
        rmidtheta(jj) = (Binstheta(jj)+Binstheta(jj+1))/2;     % Declare the average distance between bins as the representative displacement for a given bin ( centering )
    end
    
    % Save data to structure
    GLAD_binned.tsf = tsfrtheta;
    GLAD_binned.lsf = lsfrtheta;
    GLAD_binned.rdr = rmidr;  
    GLAD_binned.rdtheta = rmidtheta;
    
    GLAD_binned.counttsf = counttsfrtheta;
    GLAD_binned.countlsf = countlsfrtheta;
% 
end
