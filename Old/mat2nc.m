close all; clear all; clc;


    datapath = '/Users/JennaPalmer/Google Drive/Baylor & Jenna/Gulf of Maine/Data/Merged/';
    filename = 'test.mat';
    
    %load mat file
%     data = load([datapath filename]);
    x = ones(3,5);
    y = ones(3,5);
    
    outputdata = [datapath filename(1:end-4) '.nc'];

    nccreate(outputdata,'x',...
          'Dimensions',{'depth',3,'time',5},...
          'Format','classic')
    ncwrite(outputdata,'x',x)
    ncwriteatt(outputdata, 'x', 'long_name', 'x-dir')
    ncwriteatt(outputdata, 'x', 'Units','m')
    ncwriteatt(outputdata,'/','creation_time',datestr(now))
    ncdisp(outputdata)

    ncread(outputdata,'x')


%end