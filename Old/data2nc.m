close all; clear all; clc;

datapath = 'some/path/to/data';
filename = 'test.mat';

%load mat file
%     data = load([datapath filename]); % comment this in if you have
%     actual data to try out


% generate some fake data if you don't have a file
x = ones(3,5);
y = ones(3,5);

%outputdata = [datapath filename(1:end-4) '.nc'];
outputdata = [filename(1:end-4) '.nc'];

nccreate(outputdata,'x',...
      'Dimensions',{'depth',3,'time',5},...
      'Datatype',class(x),...
      'FillValue',NaN,...
      'Format','classic')
ncwrite(outputdata,'x',x)
ncwriteatt(outputdata, 'x', 'long_name', 'x-coor')
ncwriteatt(outputdata, 'x', 'Units','m')
ncwriteatt(outputdata,'/','creation_time',datestr(now))

nccreate(outputdata,'y',...
      'Dimensions',{'depth',3,'time',5},...
      'Datatype',class(y),...
      'FillValue',NaN,...
      'Format','classic')
ncwrite(outputdata,'y',y)
ncwriteatt(outputdata, 'y', 'long_name', 'y-coor')
ncwriteatt(outputdata, 'y', 'Units','m')
ncwriteatt(outputdata,'/','creation_time',datestr(now))
ncdisp(outputdata)

ncread(outputdata,'x')
ncread(outputdata,'y')
