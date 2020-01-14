function frame2gif(filename,cf,fn)
% Written by Jenna Pearson 2019
% filename = to include path if desired
% cf = gcf
% fn = frame number


% get frame
frame = getframe(cf); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 

% Write to the GIF File 
      if fn == 1
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 

end

