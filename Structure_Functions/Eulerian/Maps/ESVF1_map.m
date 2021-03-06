function struc = ESVF1_map(dul,dut,lon,lat,rdist,oceantime,par,outfn)
% get sizes
[~,nt] = size(lon);

if  ~isfield(par,'latlonlimbox') 
    % get lat lon bounds for bins
    maxlat = nanmax(lat(:));
    minlat = nanmin(lat(:));
    maxlon = nanmax(lon(:));
    minlon = nanmin(lon(:));
else
    % get user input lat lon bounds for bins
    maxlat = par.latlonlimbox(1,2);
    minlat = par.latlonlimbox(1,1);
    maxlon = par.latlonlimbox(2,2);
    minlon = par.latlonlimbox(2,1);
end

% Define the bin width
if  isfield(par,'binwid') 
    binwid = par.binwid;
else 
    binwid = .5; % default is .5 degree (~ 55km)
end


%-------------------------------------------------------------------------
%   Bin Vector
%-------------------------------------------------------------------------
lonbins = minlon:binwid:maxlon;
latbins = minlat:binwid:maxlat;

%-------------------------------------------------------------------------
%   Initialize matrices
%-------------------------------------------------------------------------
ul = cell([length(latbins),length(lonbins)]);
ul_over_time = cell([length(latbins),length(lonbins),length(oceantime)]);
ut = cell([length(latbins),length(lonbins)]);
ut_over_time = cell([length(latbins),length(lonbins),length(oceantime)]);
bincounts = zeros([length(latbins),length(lonbins)]);
bincounts_over_time = zeros([length(latbins),length(lonbins),length(oceantime)]);
r_over_time = cell([length(latbins),length(lonbins),length(oceantime)]);
r = cell([length(latbins),length(lonbins)]);
%-------------------------------------------------------------------------
%   Bin Data
%-------------------------------------------------------------------------
for tt = 1:nt
    for xx = 1:length(lonbins)
        for yy = 1:length(latbins)
            lontmp = lon(:,tt);
            lattmp = lat(:,tt); 
            ultmp = dul(:,tt);
            uttmp = dut(:,tt);
            rtmp = rdist(:,tt);
            
            if xx == 1 && yy == 1
                ind = find(abs(lontmp - lonbins(xx)) <= abs(lontmp - lonbins(xx+1)) & abs(lattmp - latbins(yy)) <= abs(lattmp - latbins(yy+1)));
                
                
            elseif xx == length(lonbins) && yy == length(latbins)
                ind = find(abs(lontmp - lonbins(xx)) <= abs(lontmp - lonbins(xx-1)) & abs(lattmp - latbins(yy)) <= abs(lattmp - latbins(yy-1)));
                
            elseif xx == 1 && yy == length(latbins)
                ind = find(abs(lontmp - lonbins(xx)) <= abs(lontmp - lonbins(xx+1))...
                   & abs(lattmp - latbins(yy)) <= abs(lattmp - latbins(yy-1)));
               
              
            elseif yy == 1 && xx == length(lonbins)
                ind = find(abs(lontmp - lonbins(xx)) <= abs(lontmp - lonbins(xx-1))...
                    & abs(lattmp - latbins(yy)) <= abs(lattmp - latbins(yy+1)));
                
            elseif xx == 1
                ind = find(abs(lontmp - lonbins(xx)) <= abs(lontmp - lonbins(xx+1))...
                    & abs(lattmp - latbins(yy)) <= abs(lattmp - latbins(yy+1))...
                    & abs(lattmp - latbins(yy)) <= abs(lattmp - latbins(yy-1))); 
            elseif yy == 1
                ind = find(abs(lontmp - lonbins(xx)) <= abs(lontmp - lonbins(xx+1))...
                    & abs(lontmp - lonbins(xx)) <= abs(lontmp - lonbins(xx-1))...
                    & abs(lattmp - latbins(yy)) <= abs(lattmp - latbins(yy+1)));
                
                
            elseif xx == length(lonbins)
                ind = find(abs(lontmp - lonbins(xx)) < abs(lontmp - lonbins(xx-1))...
                    & abs(lattmp - latbins(yy)) <= abs(lattmp - latbins(yy+1))...
                    & abs(lattmp - latbins(yy)) <= abs(lattmp - latbins(yy-1)));
            elseif yy == length(latbins)
                ind = find(abs(lontmp - lonbins(xx)) <= abs(lontmp - lonbins(xx+1))...
                    & abs(lontmp - lonbins(xx)) <= abs(lontmp - lonbins(xx-1))...
                    & abs(lattmp - latbins(yy)) <= abs(lattmp - latbins(yy-1)));
                
                
            else
                ind = find(abs(lontmp - lonbins(xx)) <= abs(lontmp - lonbins(xx+1))...
                    & abs(lontmp - lonbins(xx)) <= abs(lontmp - lonbins(xx-1))...
                    & abs(lattmp - latbins(yy)) <= abs(lattmp - latbins(yy+1))...
                    & abs(lattmp - latbins(yy)) <= abs(lattmp - latbins(yy-1)));   
            end

            
            if ~isempty(ind)
                if isempty(ul{yy,xx})
                    ul{yy,xx} = ultmp(ind);
                    ul_over_time{yy,xx,tt} = ultmp(ind);
                    ut{yy,xx} = uttmp(ind);
                    ut_over_time{yy,xx,tt} = uttmp(ind);
                    r{yy,xx} = rtmp(ind);
                    r_over_time{yy,xx,tt} = rtmp(ind);
                    bincounts(yy,xx) = length(ind);
                    bincounts_over_time(yy,xx,tt) = length(ind);
                else
                    ul{yy,xx}= [ul{yy,xx}; ultmp(ind)];
                    ul_over_time{yy,xx,tt} = [ul_over_time{yy,xx,tt}; ultmp(ind)];
                    ut{yy,xx} = [ut{yy,xx}; uttmp(ind)];
                    ut_over_time{yy,xx,tt} = [ut_over_time{yy,xx,tt}; uttmp(ind)];
                    r{yy,xx} = [r{yy,xx}; rtmp(ind)];
                    r_over_time{yy,xx,tt} = [r_over_time{yy,xx,tt}; rtmp(ind)];
                    bincounts(yy,xx) = length(ul{yy,xx});
                    bincounts_over_time(yy,xx,tt) = length(ul_over_time{yy,xx,tt});
                end

            end
        end
    end
    
    prev=round((tt-1)/nt,2);
    if round(tt/nt,2)~=prev
        body(['PROGRESS: ', int2str(round(tt/nt,2)*100) '%%'])
    end
end

%-------------------------------------------------------------------------
%   Statistics 
%-------------------------------------------------------------------------
ul_ave = cellfun(@nanmean,ul);
ul_sd = cellfun(@nanstd,ul);
ut_ave = cellfun(@nanmean,ut);
ut_sd = cellfun(@nanstd,ut);
r_ave = cellfun(@nanmean,r);
r_sd = cellfun(@nanstd,r);

%-------------------------------------------------------------------------
%   Boostrapping 
%-------------------------------------------------------------------------
% count = 0;
% total = length(latbins)*length(lonbins);
%  for yy = 1:length(latbins)   
%     for xx = 1:length(lonbins)
%         tmp = ul{yy,xx};
%         % make sure there are enough drifters to get a good ci and they aren't all nans
%         if length(tmp) > 8 && length(tmp(~isnan(tmp))) > 8
%             [cil,] = bootci(par.bootsampno,@nanmean,ul{yy,xx},'UseParallel','True');    
%             ulci(yy,xx,:) = cil;                       
% 
%             [cit,] = bootci(par.bootsampno,@nanmean,ut{yy,xx},'UseParallel','True');
%             utci(yy,xx,:) = cit ; 
% 
%             [cir,] = bootci(par.bootsampno,@nanmean,r{yy,xx},'UseParallel','True');
%             rci(yy,xx,:) = cir ; 
%             
%             prev=round((count-1)/total,2);
%             if round((count)/total,2)~=prev
%                 body(['PROGRESS: ', int2str(round((count)/total,2)*100) '%%'])
%             end
%             count = count+1;
%         end
%     end
%  end


%-------------------------------------------------------------------------
%   Save Data to Structure
%-------------------------------------------------------------------------

% save structures to individual files if outfn is provided including path
if nargin > 7
    struc.ul_over_time = ul_over_time;
    struc.ut_over_time = ut_over_time;
    struc.bincounts_over_time = bincounts_over_time;
    struc.r_over_time = r_over_time;
    struc.lonbins = lonbins;
    struc.latbins = latbins;
    struc.binwid = binwid;
    struc.oceantime = oceantime;

    outfn_new = [outfn(1:end-4)  '_over_time.mat'];
    save(outfn_new,'-struct','struc','-v7.3');
    clear struc

    struc.ul = ul;
    struc.ut = ut;
    struc.bincounts = bincounts;
    struc.r = r;
    struc.lonbins = lonbins;
    struc.latbins = latbins;
    struc.binwid = binwid;
    struc.oceantime = oceantime;
    
    outfn_new = [outfn(1:end-4)  '_combined.mat'];
    save(outfn_new,'-struct','struc','-v7.3');
    clear struc
    
    struc.ul_ave = ul_ave;
    struc.ul_sd = ul_sd;
%    struc.ulci = ulci;
    struc.ut_ave = ut_ave;
    struc.ut_sd = ut_sd;
%    struc.utci = utci;
    struc.r_ave = r_ave;
    struc.r_sd = r_sd;
%    struc.rci = rci;
    struc.lonbins = lonbins;
    struc.latbins = latbins;
    struc.binwid = binwid;
    struc.oceantime = oceantime;
    
    outfn_new = [outfn(1:end-4)  '_stats.mat'];
    save(outfn_new,'-struct','struc','-v7.3');
    clear struc 
end

% return the structure with all variables to manipulate without reloading
% since they are very large
struc.ul_over_time = ul_over_time;
struc.ul = ul;
struc.ul_ave = ul_ave;
struc.ul_sd = ul_sd;
%struc.ulci = ulci;

struc.ut_over_time = ut_over_time;
struc.ut = ut;
struc.ut_ave = ut_ave;
struc.ut_sd = ut_sd;
%struc.utci = utci;

struc.bincounts = bincounts;
struc.bincounts_over_time = bincounts_over_time;

struc.lonbins = lonbins;
struc.latbins = latbins;

struc.r_over_time = r_over_time;
struc.r = r;
struc.r_ave = r_ave;
struc.r_sd = r_sd;
%struc.rci = rci;

struc.binwid = binwid;

struc.oceantime = oceantime;

end