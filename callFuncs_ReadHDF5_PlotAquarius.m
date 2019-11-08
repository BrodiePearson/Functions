% ***********  Call Script and Example Usage of generic HDF5 reader and Aquarius plot Functions  ***********
%  
%  SCRIPT: callFuncs_ReadHDF5_PlotAquarius.m
%  
%  DESCRIPTION
%   - Illustrates how the HDF5 reader function and Aquarius-specific plotting functions are called from within a simple Matlab script
%   - This script can be used to read any arbitrary HDF5 file (including Aquarius L2 and L3 data files), display its metadata,
%       optionally list all file meta/data elements,and plot specified Aquarius L2 swath or L3 gridded variables.
%
%  USAGE
%   Run the script either from the Matlab command prompt or by executing via the GUI.  The script minimally invokes the 
%   data reader function to read all HDF5 file meta/data to memory, and optionally provide listings and/or call the 
%   plotting function to graphically display the data according to user selections.
%
%  DEPENDENCIES
%   The script calls the following source files and functions written in Matlab and provided by PO.DAAC:
%       readHDF5.m          readHDF5(fpath,fname,ListFlag)
%       plotAquarius.m      plotAquarius(VarData, LevelFlag, PlotParameter, TitleText, PlotMinValue, PlotMaxValue)
%   Note: see these source .m for complete function documentation.
%
%  INPUTS   ("User Inputs" code block section)
%       fpath: file path
%       fname: file name of source HDF5 data file
%       ListFlag: True/False setting to optionally output complete metadata & data listings sequentially for each file element
%       PlotFlag: True/False setting to plot Aquarius L2 or L3 data per user defined value range
%       PlotMinValue, PlotMaxValue: values defining the preferred value range for plotting
%       L2PlotParameter: (L2 only) argument from the list below (based on the dataset index within the HDF file) specifying the L2 variable to plot
%           --- Aquarius L2 dataset indices for plotting (from Group1 '/Aquarius Data') ---
%               SSS (Sea Surface Salinity)      = 22
%               Scatterometer Wind Speed        = 170
%               ancSSS (ARGO Ancillary SSS)     = 42
%               ancSST (Ancillary SST)          = 50    (in Kelvin: adjust PlotMin/Max range constants accordingly)
%               ancWindSpeed (NCEP Wind Speed)  = 57
%            
%  OUTPUTS
%       numGlobalAttributes: variable with a count of the total number of global file Attributes
%       numDatasets: variable with a count of the total number of file Datasets
%       numGroups: variable with a count of the total number of file Groups
%       globalMetadata: structure containing global file attributes
%       dataset: structure variable (dimensioned by numDatasets) of the form dataset(n).Name|Attributes|Data, where
%           Name: element name
%           Attributes: concatenated element attributes (types & values)
%           Data: element data array
%       group: structure variable (dimensioned by numGroups) of the form 
%           group(m).Name|numAttributes|Attributes|numDatasets
%           group(m).dataset(n).Name|Attributes|Data, where ...
%           numAttributes: number of element attributes
%           numDatasets: number of element datasets
%           Name: element name
%           Attributes: concatenated element attributes (types & values)
%           Data: element data array
%       [Optional] - complete listing of all HDF5 file metadata and data elements to screen
%       [Optional] - plot figure to screen
%  
%  NOTES
%      1. This read software was created using Matlab version 7.14
%      2. Please email all comments and questions concerning these routines
%           to podaac@podaac.jpl.nasa.gov.
%
%  CREATED:
%       7/11/2013: Vardis Tsontos, PO.DAAC, NASA-JPL CalTech
%
%  UPDATED:
%       10/30/2017: Vardis Tsontos, PO.DAAC, NASA-JPL CalTech
%
%  VERSION & REVISION:
%       v1.0    (001)
%
%======================================================================
% Copyright (c) 2013, California Institute of Technology
%======================================================================



% ****** User Inputs Section ***********            % User inputs entered here
Path = 'C:\Users\vtsontos\Desktop\';                % Enter the data directory path here (for both L3 and L2); inlude terminal "\" in the path string
FileName = 'Q2011362011100.L2_SCI_V5.0';            % Enter the data filename to be read (for both L3 and L2) 
ListFlag = false;                                   % Output data listings: specify true or false (for both L3 and L2)
PlotFlag = false;                                    % Output data plot: specify true or false
PlotMinValue = 30;                                  % Specify minimum value for plotting
PlotMaxValue = 40;                                  % Specify maximum value for plotting
L2PlotParameter = 22;                               % Specify L2 variable index for plotting (see INPUTS section above)
% ---------------------------------------------------------------------------------------------------------------------


% Calls the 'readAquarius()' function and returns all HDF5 file Metadata and Data as variables/structures
[numGlobalAttributes,numDatasets,numGroups,globalMetadata,dataset,group, groupInfo] = readHDF5(Path,FileName,ListFlag);  
       
% Plots the data as necessary
if (PlotFlag == true)
    % Determines which product level and parameter & dynamically constructs the plot title
    if strfind(FileName,'L3m')
        LevelFlag = 3;
        if strfind(FileName,'SSS')
            PlotParameter = 22;  %SSS
        else
            PlotParameter = 170;  %Wind Speed
        end
        % Construct the L3 Plot Title dynamically
        RegExpr= '[Q0-9]+.';
        Temp = regexp(FileName ,RegExpr,'match');
        DateTime = char(Temp(1));
        if (PlotParameter == 10)
            TitleText = ['Aquarius L3 SSS (PSU) : '  DateTime];
        else
            TitleText = ['Aquarius L3 Wind Speed (m/s) : '  DateTime];
        end;
        clear('RegExpr','Temp','DateTime');
        % Assign the variable for passing to the plotting function
        VarData = dataset;
    else
        LevelFlag = 2;
        PlotParameter = L2PlotParameter;
        indx_Orbit = 40;    %   v2.0 L2 Aquarius Global Attribute index for Orbit# 
        indx_Units = 4;     %   v2.0 L2 Aquarius Group5 '/Navigation' Parameter units index
        
        % Construct the L2 Plot Title dynamically
        RegExpr= '[Q0-9]+.';
        Temp = regexp(FileName ,RegExpr,'match');
        DateTime = char(Temp(1));
        Orbit = [char(num2str(globalMetadata(indx_Orbit).Attribute)) '  '];
        ParamName = char(groupInfo(1,1).Datasets(PlotParameter,1).Name);        
        ParamUnits = ['  (' char(groupInfo(1,1).Datasets(PlotParameter,1).Attributes(indx_Units,1).Value)  ')  '];
        TitleText = ['Aquarius L2 ' ParamName ParamUnits Orbit DateTime]; 
        clear('indx_Orbit','indx_Units','RegExpr','Temp','DateTime','ParamName','ParamUnits','Orbit');
        % Assign the variable for passing to the plotting function
        VarData = group;
    end
    % Call the plotting function
    [plot] = plotAquarius(VarData, LevelFlag, PlotParameter, TitleText, PlotMinValue, PlotMaxValue);
    
    % Final cleanup of unessential variables
    clear('plot','TitleText','VarData','groupInfo');
end;

% *******************************    END   *******************************


  



        

 
%  
        
%--------------------------------------------------------------------------
