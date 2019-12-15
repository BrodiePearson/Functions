function struc = ESVF1_map(dul,dut,lon,lat,rdist,oceantime,par,latlonlimbox)
% get sizes
[~,nt] = size(lon);

if nargin < 8
    % get lat lon bounds for bins
    maxlat = nanmax(lat(:));
    minlat = nanmin(lat(:));
    maxlon = nanmax(lon(:));
    minlon = nanmin(lon(:));
else
    % get user input lat lon bounds for bins
    maxlat = latlonlimbox(1,2);
    minlat = latlonlimbox(1,1);
    maxlon = latlonlimbox(2,2);
    minlon = latlonlimbox(2,1);
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
                % find difference between current tags and previous tags
                if isempty(ul_over_time{yy,xx})
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
                    r{yy,xx} = [r{xx,yy}; rtmp(ind)];
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

ul_ave = cellfun(@nanmean,ul);
ut_ave = cellfun(@nanmean,ut);
r_ave = cellfun(@nanmean,r);


%-------------------------------------------------------------------------
%   Boostrapping
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
%   Save Data to Structure
%-------------------------------------------------------------------------
struc.ul_over_time = ul_over_time;
struc.ul = ul;
struc.ul_ave = ul_ave;

struc.ut_over_time = ut_over_time;
struc.ut = ut;
struc.ut_ave = ut_ave;

struc.bincounts = bincounts;
struc.bincounts_over_time = bincounts_over_time;

struc.lonbins = lonbins;
struc.latbins = latbins;

struc.r_over_time = r_over_time;
struc.r = r;
struc.r_ave = r_ave;

struc.binwid = binwid;

struc.oceantime = oceantime;
end