function moviemaker(figpath,moviepath)
%   figpath = str path with location to images
%   filename = str name of movie file to be outputted without file  extension
%   ns = if there are errors, this probably needs to be adjusted to account
%   for hidden files as well as other folders.
 
    % load the images
    datanames = dir([figpath '*.png'])  ; % get list of all filenames in the data directory
    datanames = struct2cell(datanames);
    datanames = natsortfiles(datanames(1,:));
    
    % get number of images
    ni = length(datanames);
    
    images = cell(ni,1);
    
    for ii = 1:ni
        images{ii} = imread([figpath datanames{ii}]);
     end
     
     % create the video writer with 30 fps
     writerObj = VideoWriter(moviepath);
     writerObj.FrameRate = 14;
     
       % open the video writer
       open(writerObj);
       % write the frames to the video
       for ii=1:ni    
           % convert the image to a frame
           frame = im2frame(images{ii});
           writeVideo(writerObj, frame);

       end
       % close the writer object
       close(writerObj); 
end