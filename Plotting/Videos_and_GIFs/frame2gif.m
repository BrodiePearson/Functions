function frame2gif(filename,cf,fn,par)
% Written by Jenna Pearson 2019
% filename = to include path if desired
% cf = gcf
% fn = frame number
% par = structure with the following defined fields
        % par.delay = seconds between frames

% set default seconds of delay between frames (1-655)
if ~isfield(par,'delaytime')
    par.delay = 0.01;
end

% set default quality of image (1-100)
if ~isfield(par,'quality')
    par.quality = 100;
end

% get frame
frame = getframe(cf); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 

% Write to the GIF File 
      if fn == 1
          imwrite(imind,cm,filename,'gif',  'DelayTime', par.delay,...
               'Quality',par.quality, 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif', 'DelayTime', par.delay,...
              'Quality',par.quality, 'WriteMode','append'); 
      end 

end

