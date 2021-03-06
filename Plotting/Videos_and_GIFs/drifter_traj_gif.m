function drifter_traj_gif(filename,lat,lon,oceantime,varargin)
    % Written by Jenna Pearson 2020
    
    % Plots trajectories on a map background. Requires the m_map package:
    % https://www.eoas.ubc.ca/~rich/map.html. Please add this folder to
    % your matlab path prior to running.
    
    %   Required Arguments:
    %       filename:
    %       lat:
    %       lon:
    %       oceantime:
    
    
    %   Optional Arguments:
    %       color:          color of drifter head | (default black)
    %       tail_color:     color of drifter tails | (default black)
    %       tail_len:       length of drifter tail | (default 2 weeks)
    %       size_and_pos:   [left bottom height width] | (default [45, 1000000, 1300, 1300] )
    %       latlim:         [minlat maxlat] of frame | (default defined by input lat)
    %       lonlim:         [minlon maxlon] of frame | (default defined by input lon)
    %       timestep:       skip times, e.g. 4 --> 1:4:totaltime | (default 1)
    %       delay:          delay between images | (default 0.01s)

    %**********************************************************************
    % Variable Defaults
    %**********************************************************************
    
    % setup input parser to accept name value pairs
    p = inputParser;
    

    % defaults
    
    % color of drifter
    default_color = 'k'; 
    
    % tail color
    default_tail_color = 'k';
    
    % find minutes between times.
    d = (oceantime(2)-oceantime(1))*24*60;
    tail_len = 20160/d;
    
    % figure size and position
    size_and_pos = [45, 1000000, 1300, 1300];
    
    % latlim bounds of frame
    latlim = [min(lat(:)) max(lat(:))]; % define constant lat lon limits for fig
    lonlim = [min(lon(:)) max(lon(:))];
    
    % timestep
    timestep = 1;
    
    % delay between frames
    delay = 0.01;
    
    
    
    % add them to input parser
    errmsg = ['Color must be a predfined matlab character color, e.g. "k" or',...
    a 'vector of length 3.'];
    valid_func = @(x) assert(ischar(x) || isnumeric(x) && length(x) == 3,errorMsg);
    addParameter(p,'color',default_color,valid_func);
    
    errmsg = ['Color must be a predfined matlab character color, e.g. "k" or',...
    a 'vector of length 3.'];
    valid_func = @(x) assert(ischar(x) || isnumeric(x) && length(x) == 3,errorMsg);
    addParameter(p,'tail_color',default_tail_color,valid_func);
    
    parse(p,varargin{:});
    
    font = p.Results.font;
    fontsize = p.Results.fontsize;
    
    % set values in p to variables
    
%     if ~exist('color','var')
%         color = 'k';
%     end
    
   % color of drifter tail
    
%     if ~exist('tail_color','var')
%         color = 'k';
%     end
    
    % length of drifter tail
%     if ~exist('tail_len','var')
%     end
%     
%     % figure size and position
%     if ~exist('size_and_pos','var')
%         size_and_pos = [45, 1000000, 1300, 1300];
%     end
%     
%     % latlim bounds of frame
%     if ~exist('latlim','var') && ~exist('lonlim','var')
%         latlim = [min(lat(:)) max(lat(:))]; % define constant lat lon limits for fig
%         lonlim = [min(lon(:)) max(lon(:))];
%     end
%     
%     % timestep
%     if ~exist('timestep','var')
%         timestep = 1;
%     end
%     
%     % delay between frames
%     if ~exist('delay','var')
%         delay = 0.01;
%     end
%     

    %**********************************************************************
    % Main Code
    %**********************************************************************
    
    figure;
    set(gcf,'position',size_and_pos) %  set figure position and size
    axis tight manual

    for tt = 1:timestep:length(oceantime)

        m_proj('lambert','long',lonlim,'lat',latlim,'rectbox','on');  
        
        if tt < tail_len
            tmplat = lat(:,1:tt-1); tmplon = lon(:,1:tt-1);
            m_scatter(tmplon(:),tmplat(:),800,'k.')
        else
            m_scatter(lon(:,tt-tail_len:tt-1),lat(:,tt-tail_len:tt-1),800,tail_color,'.')
        end
        
        hold on

        m_scatter(lon(:,tt),lat(:,tt),1000,color,'.')

        m_grid('box','on','tickdir','in','backgroundcolor',clr('light_blue')./255);
        m_gshhs_i('color','k');  
        m_gshhs_i('patch',[.7 .7 .7]);


        title(['Trajectories '  datestr(oceantime(tt))])
        xlabel('Longitude')
        ylabel('Latitude')

        pause(.5)
        frame2gif(filename,gcf,tt)
        hold off
    end


end