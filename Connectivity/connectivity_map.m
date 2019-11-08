function struc = connectivity_map(lon,lat,oceantime,par)
% get sizes
[nd,nt] = size(lon);

% creat tags for each drifter
tag = 1:nd;

% get lat lon bounds for boxes
maxlat = max(lat(:));
minlat = min(lat(:));
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
first_hit = cell([length(latbins),length(lonbins)]);
first_hit_time = cell([length(latbins),length(lonbins)]);
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
                cur_tags = tag(ind);
                cur_time = oceantime(tt);
                % find difference between current tags and previous tags
                if isempty(first_hit{yy,xx})
                    first_hit{yy,xx} = cur_tags;
                    first_hit_time{yy,xx} = cur_time.*ones(size(cur_tags));
                    bincounts(yy,xx) = length(ind);
                    bincounts_over_time(yy,xx,tt) = length(ind);
                else
                    old_tags_ind = ismember(tag(ind),first_hit{yy,xx});
                    new_tags = cur_tags(~old_tags_ind);
                    first_hit{yy,xx} = [first_hit{yy,xx} new_tags];
                    first_hit_time{yy,xx} = [first_hit_time{yy,xx} cur_time.*ones(size(new_tags))];
                    bincounts(yy,xx) = length(first_hit{yy,xx});
                    bincounts_over_time(yy,xx,tt) = length(first_hit{yy,xx});
                end

            end
        end
    end
    
    prev=round((tt-1)/nt,2);
    if round(tt/nt,2)~=prev
        body(['PROGRESS: ', int2str(round(tt/nt,2)*100) '%%'])
    end
end

total_prob = bincounts./nd;

%-------------------------------------------------------------------------
%   Save Data to Structure
%-------------------------------------------------------------------------
struc.bincounts = bincounts;
struc.bincounts_over_time = bincounts_over_time;
struc.total_prob = total_prob;
struc.lonbins = lonbins;
struc.latbins = latbins;
struc.binwid = binwid;
struc.oceantime = oceantime;
struc.first_hit = first_hit;
struc.first_hit_time = first_hit_time;
end