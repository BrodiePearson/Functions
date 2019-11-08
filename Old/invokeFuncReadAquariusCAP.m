% ***********  Call Script and Example Usage of Aquarius-CAP L2 & L3 Read Functions  ***********
%  
%  SCRIPT: invokeFuncReadAquariusCAP.m
%  
%  DESCRIPTION
%   - Illustrates how Aquarius-CAP L2 and L3 functions are called from within a simple Matlab script
%   - This script can be used to read Aquarius-CAP L2 and L3 data files, display metadata & data listings,
%       and create plots of specified variables.
%
%  USAGE
%   The script automatically invokes either of the following Aquarius-CAP L2 or L3 reader functions depending on the 
%       file-type the user specifies as input
%
%  DEPENDENCIES
%   The script invokes the following CAP reader functions and source files written in Matlab and provided by PO.DAAC:
%       readL2AquariusCAP.m   readL3AquariusCAP()
%       readL3AquariusCAP.m   readL2AquariusCAP()
%   Note: see the source L2 and L3 .m files for complete function documentation.
%
%  INPUTS   ("User Inputs" code block section)
%       fpath: file path   (required for both L3 & L2)
%       fname: file name of source Aquarius-CAP L2 HDF5 or L3 NetCDF data file    (required for both L3 & L2)
%       DisplayListings: True/False setting to output data listings sequentially by Dataset parameter (required for both L3 & L2)
%       DisplayPlot: True/False setting to plot mapped L3 SSS data within a user defined range (required only for L3)
%       minSalinityPlotValue, maxSalinityPlotValue: values defining range of salinity values plotted (required only for L3)
%       indxPlotDataset: arguement from the list belowspecifying whether L2 data is plotted and for which variable (L2 only).
%       CAP-L2 plot indices and data variables:
%           1= 'SSS'                  2= 'SSS_CAP'             3= 'SSS_CAP_h'
%           4= 'SSS_CAP_v'            5= 'anc_SSS'             6= 'anc_surface_temp'
%           7= 'anc_wind_dir'         8= 'anc_wind_speed'      9= 'beam_clat'
%           10= 'beam_clon'           11= 'cap_flag'           12= 'ice_frac'
%           13= 'land_frac'           14= 'scat_ice_frac'      15= 'scat_land_frac'
%           16= 'scat_wind_speed'     17= 'sec'                18= 'wind_dir_cap'
%           19= 'wind_speed_cap'       
%           0= No plot will be generated (equivalent to DisplayFlag = False)
%            
%  OUTPUTS
%       globalDimensions: structure containing global dimension attributes for arrays
%       globalMetadata: structure containing global file attributes
%       VarData: multidimensional data structure (dimensioned on #variables)
%       VectorData: (L3 only) CAP-Salinity or WindSpeed data in vector array (Lat, Lon, SSS) format
%       optional plot figure to screen
%       optional Listings of dataset attributes & data to screen
%  
%  NOTES
%      1. This read software was created using Matlab version 7.14
%      2. Please email all comments and questions concerning these routines
%           to podaac@podaac.jpl.nasa.gov.
%
%  CREATED:
%       7/26/2012: Vardis Tsontos, PO.DAAC, NASA-JPL CalTech
%
%  UPDATED:
%       4/09/2013: Vardis Tsontos, PO.DAAC, NASA-JPL CalTech
%       10/21/2013:  Vardis Tsontos, PO.DAAC, NASA-JPL CalTech  (generalize for CAP-CF version)       
%       
%
%======================================================================
% Copyright (c) 2013, California Institute of Technology
%======================================================================



% ****** User Inputs Section ***********            % User inputs entered here
Path = 'C:\Users\vtsontos\Desktop\';                % Enter the data directory path here (for both L3 and L2)
Name = 'Q2011335193600.L2_SCI_V3.0.cap';            % Enter the data filename to be read (for both L3 and L2) 
ListFlag = false;                                   % Output data listings: specify true or false (for both L3 and L2)
DisplayFlag = true;                                 % Output data plot: specify true or false (L3 only)
minValue = 30;                                      % Specify minimum value for plotting (L3 only)
maxValue = 40;                                      % Specify maximum value for plotting (L3 only)
indxPlotDataset = 2;                                % Specify L2 variable index for plotting (see INPUTS section above)
% ---------------------------------------------------------------------------------------------------------------------



% *******************************    Main   *******************************

% Determines which product type to process & underlying function to invoke based on specified filename pattern (Default:L2)
if strfind(Name,'cap.nc')
    LevelFlag = 3;
    if strfind(Name,'sss')
        L3parameter = 'Salinity';
    else
        L3parameter = 'Wind Speed';
    end
else
    LevelFlag = 2;
end;

if LevelFlag == 3
    [globalDimensions, globalMetadata, VarData VectorData]=readL3AquariusCAP(Path,Name,ListFlag,DisplayFlag,L3parameter, minValue,maxValue); %L3 function call    
else
    [globalMetadata, VarData] = readL2AquariusCAP(Path, Name, ListFlag, indxPlotDataset);              %L2 function call
end;

%--------------------------------------------------------------------------
