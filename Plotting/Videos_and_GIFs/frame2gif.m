function frame2gif(filename,cf,fn,varargin)
% Written by Jenna Pearson 2019
% filename = to include path if desired
% cf = gcf
% fn = frame number

% Optional Arguments
%       delay = seconds between frames }| (default 0.01s)

% set default seconds of delay between frames (1-655)
if ~exist('delay','var')
    delay = 0.01;
end


% get frame
frame = getframe(cf); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 

% Write to the GIF File 
      if fn == 1
          imwrite(imind,cm,filename,'gif',  'DelayTime', delay,...
                'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif', 'DelayTime', delay,...
               'WriteMode','append'); 
      end 

end

