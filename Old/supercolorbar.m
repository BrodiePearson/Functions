function supercolorbar(gcf,orientation,label,fontsize, min, max)
% creates a single colorbar with uniform max and mins and applies only to the final subplot
% in a series keeping all the figure positions and dimensions intact.

% Inputs:
%   label = label for caxis ex: 'm/s'
%   location = string value of placement and orientation  ex: 'east' or 'south'
%   min = minimum value for colorbar;
%   max = maximum value for colorbar;

% If no min and max supplied it defaults to the largest max and smallest min
% close all; clear all; clc;

% figure(1)
% subplot(1,2,1)
% quiver(1:5,1:5)
% 
% subplot(1,2,2)
% m2 = magic(5).*.7;
% pcolor(m2)
%     
% subplot(1,3,1)
% m1 = magic(5);
% pcolor(m1)
% subplot(1,3,2)
% m2 = magic(5).*.7;
% pcolor(m2)
% subplot(1,3,3)
% m2 = magic(5).*.7;
% pcolor(m2)

% orientation = 'south';
% label = '(m/s)^2';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orientation = [orientation 'outside'];

% grab axes for all subplots, note the first one is the last subplot
h=get(gcf,'children');
gca = get(h,'type');


% defaults when some arguments are not submitted.
if nargin < 5
    data=findobj(gcf,'type','axes');
    
    for dd = 1:length(data)
        f=get(h(dd),'children');
        t = get(f,'type');
        
        if isequal(t,'surface')
            ctemp = get(f,'cdata');
            ctemp = ctemp(:);
            cmin(dd)=nanmin(ctemp);
            cmax(dd)=nanmax(ctemp);
        end
    end
    
    miny = nanmin(cmin);
    maxy = nanmax(cmax);
elseif nargin < 4
    fontsize = 20;
end
    

for dd = 1:length(h)
    ax = h(dd,1);
    c = colorbar(ax,orientation);
    caxis(ax,[miny,maxy])
    c.Title.String = label;
    
    if dd > 1
        colorbar(ax,'off')
    else
        if strcmp(orientation,'eastoutside') == 1 
            p1 = get(c,'Position');
            c.Position = p1 + [.05 0 0 0]; % [left bottom width height]
        else
            p1 = get(c,'Position');
            c.Position = p1 + [-p1(1)+.25 -.1 -p1(3)+ .5 0]; % [left bottom width height]
        end
    end
end
end