function [globalMetadata, VarData] = readL2AquariusCAP(fPath, fName, DisplayListings, DisplayPlotVar)
%  
%  FUNCTION: readL2AquariusCAP()
%
%  DESCRIPTION
%   This function reads a user defined Aquarius-CAP HDF5 file and returns a data structure array (VarData) 
%   containing complete dataset attributes & data.  Input arguements require file path and name specification, and provide
%   the user with options to plot the L2 swath data and/or dump both
%   dataset attributes & data as listings to screen for review.
%
%  USAGE
%   Call the "readL2AquariusCAP() " function from either the Matlab command
%   line or within a script using suitable arguement values
%
%  INPUT ARGUMENTS
%       fpath: file path  (embed in single quotes eg.  'C:\')
%       fname: file name of source HDF5 data file (embed in single quotes eg. 'Q2013001004300.L2_SCI_V2.0.cap.hdf' )
%       DisplayListings: 0/1 (True/False) flag setting to output data listings sequentially by Dataset parameter
%       DisplayPlotVar: plot flag/variable selection: for v3.0 of the CAP
%       data (note that some indices have changed relative to v2.0 with the
%       inclusion of additional variables.  Setup here is for v3.0 but see line 129 below to adjust for v2.0 plotting)
%           0 (or any value outside range 1-9 that corresponds to a valid Dataset parameter) = plot flag is False
%           1= 'SSS'                  2= 'SSS_CAP'             3= 'SSS_CAP_h'
%           4= 'SSS_CAP_rc'           5= 'SSS_CAP_v'           6= 'anc_SSS'
%           7= 'anc_surface_temp'     8= 'anc_swh'             9= 'anc_wind_dir'
%           10= 'anc_wind_speed'      11= 'beam_clat'          12= 'beam_clon'
%           13= 'cap_far_rot'         14= 'cap_flag'           15= 'ice_frac'
%           16= 'land_frac'           17= 'scat_ice_frac'      18= 'scat_land_frac'
%           19= 'scat_swh_wind_speed' 20= 'scat_wind_speed'    21= 'sec' 
%           22= 'wind_dir_cap'        23= 'wind_speed_cap'       
%
%  OUTPUTS
%       globalMetadata: structure (1 element) containing global file attributes
%       VarData: multidimensional data structure (dimensioned on #datasets = 19)
%           Name: parameter/variable name
%           Attributes: concatenated parameter attributes (types & values)
%           Data: 3 x 4084 array with parameter/variable data for each beam
%       optional plot figure to screen
%       optional Listings of dataset attributes & data to screen
%  
%  NOTES
%      1. This read software was created using Matlab version 7.14
%      2. Please email all comments and questions concerning these routines
%           to podaac@podaac.jpl.nasa.gov.
%
%  CREATED:
%       7/23/2012: Vardis Tsontos, PO.DAAC, NASA-JPL CalTech
%
%  UPDATED:
%       4/09/2013: Vardis Tsontos, PO.DAAC, NASA-JPL CalTech
%       10/22/2014: Vardis Tsontos, PO.DAAC, NASA-JPL CalTech
%
%======================================================================
% Copyright (c) 2012, California Institute of Technology
%======================================================================



% *************************** Input Settings ***************************

%Version 2.0 CAP CAP Dataset parameters and indices:
% 1= 'SSS'                  2= 'SSS_CAP'             3= 'SSS_CAP_h'
% 4= 'SSS_CAP_v'            5= 'anc_SSS'             6= 'anc_surface_temp'
% 7= 'anc_wind_dir'         8= 'anc_wind_speed'      9= 'beam_clat'
%10= 'beam_clon'           11= 'cap_flag'           12= 'ice_frac'
%13= 'land_frac'           14= 'scat_ice_frac'      15= 'scat_land_frac'
%16= 'scat_wind_speed'     17= 'sec'                18= 'wind_dir_cap'
%19= 'wind_speed_cap' 

%Version 3.0 CAP CAP Dataset parameters and indices:
%1= 'SSS'                  2= 'SSS_CAP'             3= 'SSS_CAP_h'
%4= 'SSS_CAP_rc'           5= 'SSS_CAP_v'           6= 'anc_SSS'
%7= 'anc_surface_temp'     8= 'anc_swh'             9= 'anc_wind_dir'
%10= 'anc_wind_speed'      11= 'beam_clat'          12= 'beam_clon'
%13= 'cap_far_rot'         14= 'cap_flag'           15= 'ice_frac'
%16= 'land_frac'           17= 'scat_ice_frac'      18= 'scat_land_frac'
%19= 'scat_swh_wind_speed' 20= 'scat_wind_speed'    21= 'sec' 
%22= 'wind_dir_cap'        23= 'wind_speed_cap'


DisplayPlot = false;
if ((DisplayPlotVar > 0) && (DisplayPlotVar < 20))
    DisplayPlot = true;
    PlotParameter = DisplayPlotVar;
end


% ***************************** Main **************************************

% Global Constants
numBeams = 3;

%Read HDF5 CAP data and variable attributes into VarData structure
FileName = [fPath fName];
FileInfo = h5info(FileName);
numGlobalAttributes = length(FileInfo.Attributes);
numGroups = length(FileInfo.Groups);
numDatasets = length(FileInfo.Datasets);

globalAttributes = FileInfo.Attributes;
for n =1: numGlobalAttributes
    globalMetadata(n).Attribute = [globalAttributes(n).Name  ': ' globalAttributes(n).Value];
end

for n = 1:numDatasets
    numAttributes = length(FileInfo.Datasets(n).Attributes);
    Temp ='';
    for m = 1:numAttributes
        Temp = [Temp '   ' FileInfo.Datasets(n).Attributes(m).Name  ': ' num2str(FileInfo.Datasets(n).Attributes(m).Value)];
    end
    VarData(n).Name = FileInfo.Datasets(n).Name;
    VarData(n).Attributes = strtrim(Temp);
    VarData(n).Data = h5read(FileName,['/' FileInfo.Datasets(n).Name]);
end

% Plot the swath data for all beams
if DisplayPlot
    RegExpr= '[Q0-9]+.';
    Temp = regexp(fName ,RegExpr,'match');
    DateTime = char(Temp(1));
    ParamName = [char(FileInfo.Datasets(PlotParameter).Attributes(1).Value) '  ('];
    ParamUnits = [char(FileInfo.Datasets(PlotParameter).Attributes(2).Value)  ') :     '];
    TitleText = [ParamName ParamUnits DateTime];
    
    for n= 1:numBeams
        ParamData = VarData(PlotParameter).Data(n,:);
        indexParamDatainvalid = find(ParamData < 0);
        ParamData(indexParamDatainvalid) = NaN;
 %       scatter(VarData(10).Data(n,:),VarData(9).Data(n,:),5, ParamData);
 %       use the above scatter command (line 115) for v2.0 of the CAP L2 data and comment out the below
 %       or if working with CAP v3.0 comment out the below scatter command (line 121)
 %       and use the above one.  This accommodates changes in the Array
 %       index of the Lat and Lon variables within the structure between
 %       versions 2.0 and 3.0 of the CAP L2 dataset.
         scatter(VarData(12).Data(n,:),VarData(11).Data(n,:),5, ParamData);
        axis([-180 180 -90 90]);
        if n == 1
            xlabel('Longitude');
            ylabel('Latitude');
            title(TitleText, 'FontSize', 10);
            colorbar;
            hold on;
        end
    end   
end

% Display listing of file dataset attributes and values
if DisplayListings
    disp(' ');    
    disp('     Global Dataset Attributes');
    disp('------------------------------------');
    disp(['File Name: ' num2str(FileName)]);
    disp(['Number of Global Attributes: ' num2str(numGlobalAttributes)]);
    disp(['Number of Groups: ' num2str(numGroups)]);
    disp(['Number of Datasets: ' num2str(numDatasets)]);

    disp(' ');
    disp(['--------- Press Space to Continue with Global Metadata listings  --------------']);
    pause;
    disp(' ');
    for n = 1:numGlobalAttributes
        disp(['Attribute # ' num2str(n)]);
        disp([globalMetadata(n).Attribute]);   
        disp(' ');
    end
    
    disp(' ');
    disp(['--------- Press Space to Continue with Variable Attribute & Data listings  --------------']);
     pause;
    disp(' ');
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
clear('FileInfo', 'fPath', 'fName');
clear('m', 'n', 'Temp');
clear('indexParamDatainvalid', 'ParamData', 'numAttributes');
clear('RegExpr', 'DateTime', 'ParamName', 'ParamUnits', 'TitleText', 'PlotParameter');
clear('DisplayListings', 'DisplayPlot');

% ***************************** End **************************************


end

