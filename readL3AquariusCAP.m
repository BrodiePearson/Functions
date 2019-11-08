function [globalDimensions, globalMetadata, VarData, VectorData] = readL3AquariusCAP(fPath, fName, DisplayListings, DisplayPlot, L3parameter, varargin)
%  
%  FUNCTION: readL3AquariusCAP()
%
%  DESCRIPTION
%   This function reads a user defined Aquarius-CAP L3 NetCDF file and returns a data structure array (VarData) 
%   containing complete source dataset attributes & Sea Surface Salinity (SSS) or Wind Speed data. It also returns
%   a filtered version of the  data array in vectorized form (VectorData). Input arguements require file path and name
%   specification, and provide the user with options to plot the L3 mapped data and/or dump both
%   dataset attributes & data as listings to screen for review.
%
%  USAGE
%   Call the "readL3AquariusCAP() " function from either the Matlab command
%   line or within a script using suitable arguement values
%
%  INPUTS ARGUEMENTS
%       fname: file name of source Aquarius-CAP NetCDF data file
%       fpath: file path
%       DisplayListings: True/False flag setting to output data listings sequentially by Dataset parameter
%       DisplayPlot: True/False flag setting to create plot of mapped L3 SSS data within a user defined range
%       L3parameter: Level 3 file core data parameter ('Salinity' or 'Wind Speed')
%       minValue, maxValue: inputs defining the range of salinity or wind speed values plotted & contained in the VectorData structure.
%         (note: All source data present in VarData storage structures & listed on screen irrespective of applied filters)
%
%  OUTPUTS
%       globalMetadata: structure (10 elements) containing global file attributes
%       VarData: multidimensional data structure (dimensioned on #variables = 18)
%           Name: parameter/variable name
%           Attributes: concatenated parameter attributes (types & values)
%           Data: 360 x 181 array with source mapped 1 degree CAP Sea surface
%           salinity values (includes 'Null/NoFill' values = 32767)
%       VectorData: structure with lat/lon/SSS or WindSpeed variables in vectorized form
%           (65160 elements; range filers applied with NaN assignments)
%       optional plot figure to screen
%       optional Listings of dataset attributes & data to screen
%  
%  NOTES
%      1. This read software was created using Matlab version 7.14
%      2. Please email all comments and questions concerning these routines
%           to podaac@podaac.jpl.nasa.gov.
%
%  CREATED:
%       7/24/2012: Vardis Tsontos, PO.DAAC, NASA-JPL CalTech
%
%  UPDATED:
%       4/9/2013: Vardis Tsontos, PO.DAAC, NASA-JPL CalTech
%       10/21/2013:  Vardis Tsontos, PO.DAAC, NASA-JPL CalTech  (generalize for CAP-CF version) 
%
%======================================================================
% Copyright (c) 2013, California Institute of Technology
%======================================================================




% ***************************** Main **************************************

% Use default or user defined SSS/Wind Speed display range
if strcmp(L3parameter,'Salinity')
    ObsMin = 30;
    ObsMax = 40;
else
    ObsMin =  0;
    ObsMax = 30; 
end;

nVarargs = length(varargin);
if nVarargs == 2
   Val1 = cell2mat(varargin(1));
   Val2 = cell2mat(varargin(2));
   % assigns min/max values to correct variable irrespective of Arguement order user provides
   if(Val1 > Val2)
        ValMin = Val2;
        ValMax = Val1;
   else
       ValMin = Val1;
       ValMax = Val2;       
   end
   clear('Val1', 'Val2');
   % checks that arguement values within an acceptable range before adopted as min/max filters
   if (ValMin >= 0)
       ObsMin = ValMin;
   end
   if (ValMax <= ObsMax)
       ObsMax = ValMax;
   end
end

% Read file Dimensions metadata into GlobalDims structure
FileName = [fPath fName];
FileInfo = ncinfo(FileName);
FileFormat = FileInfo.Format;
globalDimensions = FileInfo.Dimensions;
numGlobalDimensions = length(FileInfo.Dimensions);
for n=1: numGlobalDimensions 
    if strfind(lower(globalDimensions(n).Name), 'lon')
        LonBins = globalDimensions(n).Length;
    elseif strfind(lower(globalDimensions(n).Name), 'lat')
            LatBins = globalDimensions(n).Length;        
    end;
    globalDimensions(n).Attribute = [globalDimensions(n).Name  ': ' num2str(globalDimensions(n).Length)];
end
    
% Read file metadata & global CAP dataset attributes into globalMetadata structure
globalAttributes = FileInfo.Attributes;
numGlobalAttributes = length(FileInfo.Attributes);
numGroups = length(FileInfo.Groups);
numDatasets = length(FileInfo.Variables);
for n =1: numGlobalAttributes
    if isnumeric(globalAttributes(n).Value)
            attribVal = num2str(globalAttributes(n).Value(:)');
        else
            attribVal = globalAttributes(n).Value;
        end    
        globalMetadata(n).Attribute = [globalAttributes(n).Name  ':  ' attribVal];
 end


% Read NetCDF CAP data and variable attributes into VarData structure
for n = 1:numDatasets
    numAttributes = length(FileInfo.Variables(n).Attributes);
    Temp ='';
    for m = numAttributes:-1:1
        Temp = [Temp '   ' FileInfo.Variables(n).Attributes(m).Name  ': ' num2str(FileInfo.Variables(n).Attributes(m).Value)];
    end
    VarData(n).Name = FileInfo.Variables(n).Name;
    VarData(n).Attributes = strtrim(Temp);
    VarData(n).Data = ncread(FileName,['/' FileInfo.Variables(n).Name]);
end


% preallocates VectorData structure array - optimization
VectorData.Lon = zeros(1, LonBins*LatBins);
VectorData.Lat = zeros(1, LonBins*LatBins);
VectorData.Obs = zeros(1, LonBins*LatBins);
cnt = 0;
% Vectorizes the CAP dataset and filters off out of bound values
for m = 1:LonBins
    Lon = VarData(1).Data(m);
    if Lon >= 180
        Lon = -360 + Lon;
    end
    for n = 1:LatBins
        Lat = VarData(2).Data(n);
        cnt = cnt +1;
        Obs = VarData(3).Data(m,n);
        if ((Obs < ObsMin) || (Obs > ObsMax))
            Obs = NaN;
        end
        VectorData.Lon(cnt) = Lon;
        VectorData.Lat(cnt) = Lat;
        VectorData.Obs(cnt) = Obs;
    end
end

% Plot the Gridded L3 CAP data
if DisplayPlot
    RegExpr= '[2010-9]+.';
    Temp = regexp(fName ,RegExpr,'match');
    DateTime = char(Temp(1));
    if strcmp(L3parameter,'Salinity')
        TitleText = ['Aquarius CAP Sea Surface Salinity (PSU) : '  DateTime];
    else
        TitleText = ['Aquarius CAP Wind Speed (m/s) : '  DateTime];
    end;
    scatter(VectorData.Lon,VectorData.Lat,12, VectorData.Obs, 'square', 'filled');
    axis([-180 180 -90 90]);
    title(TitleText, 'FontSize', 10);
    xlabel('Longitude');
    ylabel('Latitude');
    colorbar;
end

% Display listing of file dataset attributes and values
if DisplayListings
    disp(' ');
    disp(['File Name: ' FileName  '      NetCDF version: ' FileFormat]);
    disp(['Number of Global Dimensions: ' num2str(numGlobalDimensions)]);
    disp(['Number of Global Attributes: ' num2str(numGlobalAttributes)]);
    disp(['Number of Groups: ' num2str(numGroups)]);
    disp(['Number of Datasets: ' num2str(numDatasets)]);
    disp(' ');    
    disp('     Global Dataset Dimensions & Attributes');
    disp('------------------------------------');
    for n =1: numGlobalDimensions
        disp(num2str(globalDimensions(n).Attribute));
    end
    disp(' '); 
    for n =1:numGlobalAttributes
        disp(num2str(globalMetadata(n).Attribute));
    end
    disp(' ');
    disp(['--------- Press Space to Continue with Variable Attribute & Data listings  --------------']);
    pause;
    disp(' ');
    disp('     Variable Attributes');
    disp('------------------------------------');
    for n = 1:numDatasets
        disp(['Dataset # ' num2str(n)]);
        disp(['Dataset Name: ' VarData(n).Name]);
        disp(['Dataset Attributes: ' VarData(n).Attributes]);
        disp('Dataset Values: ');
        disp(VarData(n).Data);
        if n == numDatasets
            disp('------------------   End   -------------------');
        else
            disp(['--------- Press Space to Continue for Remaining ' num2str(numDatasets-n) ' Datasets  --------------']);
            pause;
        end       
    end
end
    

% Cleanup temporary variables
clear('FileInfo', 'fPath', 'fName','globalAttributes');
clear('m', 'n','Temp');
clear('Lat', 'Lon', 'LatBins', 'LonBins', 'Obs', 'cnt', 'numAttributes','attribVal');
clear('RegExpr', 'DateTime', 'TitleText');
clear('DisplayListings');
clear('DisplayPlot');

% ***************************** End **************************************
