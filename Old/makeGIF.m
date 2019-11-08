function makeGIF(gcf,framenum, filename,speed)

    frame = getframe(gcf); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    
    if framenum == 1
        imwrite(imind,cm,filename,'gif','DelayTime',speed, 'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename,'gif','DelayTime',speed,'WriteMode','append'); 
    end 
    
    
end