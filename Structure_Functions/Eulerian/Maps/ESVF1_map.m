function struc = ESVF1_map(dul,dut,lon,lat,oceantime,par)
% get sizes
[nd,nt] = size(lon);

% creat tags for each drifter
tag = 1:nd;

% get lat lon bounds for boxes
maxlat = max(lat(:));
minlat = nanmin(lat(:));
maxlon = max(lon(:));
minlon = min(lon(:));

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
ul_total= cell([length(latbins),length(lonbins)]);
ut_total = cell([length(latbins),length(lonbins)]);
bincounts = zeros([length(latbins),length(lonbins)]);
bincounts_over_time = zeros([length(latbins),length(lonbins),length(oceantime)]);
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
                if isempty(ul_total{yy,xx})
                    ul_total{yy,xx} = ultmp(ind);
                    ut_total{yy,xx} = uttmp(ind);
                    bincounts(yy,xx) = length(ind);
                    bincounts_over_time(yy,xx,tt) = length(ind);
                else
                    ul_total{yy,xx}= [ul_total{yy,xx}; ultmp(ind)];
                    ut_total{yy,xx} = [ut_total{yy,xx}; uttmp(ind)];
                    bincounts(yy,xx) = length(ul_total{yy,xx});
                    bincounts_over_time(yy,xx,tt) = length(ul_total{yy,xx});
                end

            end
        end
    end
    
    prev=round((tt-1)/nt,2);
    if round(tt/nt,2)~=prev
        body(['PROGRESS: ', int2str(round(tt/nt,2)*100) '%%'])
    end
end

avul = cellfun(@nanmean,ul_total);
avut = cellfun(@nanmean,ut_total);
%-------------------------------------------------------------------------
%   Save Data to Structure
%-------------------------------------------------------------------------
struc.ul_total = ul_total;
struc.ut_total = ut_total;
struc.avul = avul;
struc.avut = avut;
struc.bincounts = bincounts;
struc.bincounts_over_time = bincounts_over_time;
struc.lonbins = lonbins;
struc.latbins = latbins;
struc.binwid = binwid;
struc.oceantime = oceantime;
end