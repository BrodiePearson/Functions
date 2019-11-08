function [plot] = plotAquarius(VarData, LevelFlag, PlotParameter, TitleText, varargin)
%  
%  FUNCTION: plotAquarius(VarData, LevelFlag, PlotParameter, TitleText, PlotMinValue, PlotMaxValue)
%
%  DESCRIPTION
%   This function plots the data from a user selected Level 2 or Level 3 Aquarius HDF5 file and parameter of choice to screen 
%   for visualization of the swath or gridded data. This function merely plots the data and assumes that the HDF
%   file has been read in advance (see DEPENDENCIES below). Input arguments include the source data variable, the product 
%   level and parameter for display, with optional value range arguements for both filtering and scaling of data values.
%
%  USAGE    plotAquarius(VarData, LevelFlag, PlotParameter, TitleText, PlotMinValue, PlotMaxValue);
%   Call the "plotAquarius()" function from either the Matlab command line
%       or within a script using suitable arguement values
%
%  DEPENDENCIES
%   The Aquarius data for plotting must have been read into memory using
%       the readHDF5() function or analogue in advance
%
%  INPUT ARGUMENTS
%       VarData: the variable structure for plotting
%       LevelFlag: True/False flag setting to create plot of mapped L3 SSS data within a user defined range
%       PlotParameter: data parameter index to plot (eg. for L3 'Salinity' or 'Wind Speed', for L2 additional paremeters too)
%       TitleText: plot title text (dynamically constructed externally)
%       minValue, maxValue (Optional): inputs defining the range of values plotted (and contained in the local/temporary VectorData structure)       
%
%  OUTPUTS
%       plot of selected Aquarius L2 or L3 parameter
%  
%  NOTES
%      1. This read software was created using Matlab version 7.14
%      2. Please email all comments and questions concerning these routines
%           to podaac@podaac.jpl.nasa.gov.
%
%  CREATED:
%       7/11/2013: Vardis Tsontos, PO.DAAC, NASA-JPL CalTech
%
%  LAST UPDATED:
%       10/30/2017: Vardis Tsontos, PO.DAAC, NASA-JPL CalTech
%
%  VERSION & REVISION:
%       v1.0    (001)
%
%======================================================================
% Copyright (c) 2013, California Institute of Technology
%======================================================================

LonTag = 'Longitude';
LatTag = 'Latitude';

%Register available user-defined parameter value range or apply a default range
if (PlotParameter == 22)||(PlotParameter == 42)    % default SSS range
    ObsMin = 30;
    ObsMax = 40;
else
    ObsMin =  0;            % Assumes alternate PlotParameter is WindSpeed or SST with this default range
    ObsMax = 30; 
end 

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
    clear('ValMax', 'ValMin');
end

% Plot either the L2 of L3 Aquarius data
if (LevelFlag == 2)
    numBeams = 3;
    indx_beam_clat = 3;     % Group5 '/Navigation' Dataset indices in v2.0 of Aqurius dataset
    indx_beam_clon = 4;   
    
    for n= 1:numBeams
        ParamData = VarData(1).Dataset(PlotParameter).Data(n,:);
        indexParamDatainvalid = find((ParamData < ObsMin));
        ParamData(indexParamDatainvalid) = NaN;
        indexParamDatainvalid = find((ParamData > ObsMax));
        ParamData(indexParamDatainvalid) = NaN;
        LatData = VarData(5).Dataset(indx_beam_clat).Data(n,:);
        LonData = VarData(5).Dataset(indx_beam_clon).Data(n,:);            
        scatter(LonData,LatData,5, ParamData, 'square', 'filled');  % plots unfiltered block data values
        axis([-180 180 -90 90]);
        if n == 1
            xlabel(LonTag);
            ylabel(LatTag);
            title(TitleText, 'FontSize', 10);
            colorbar;
            hold on;
        end
    end   
else      % L3
    LonBins = 360;
    LatBins = 180;
        
    % preallocates VectorData structure array - optimization
    VectorData.Lon = zeros(1, LonBins*LatBins);
    VectorData.Lat = zeros(1, LonBins*LatBins);
    VectorData.Obs = zeros(1, LonBins*LatBins);
    cnt = 0;
    % Vectorizes the Aquarius L3 dataset and filters off out of bound values
    for m = 1:LonBins
        Lon = m - 180.5;   % Adjusts for Aquarius L3m data Array X=1 corrsponding to Dateline (-180 Lon)
        for n = 1:LatBins
            Lat = 90.5 - n;
            cnt = cnt +1;
            Obs = VarData(1).Data(m,n);
            if ((Obs < ObsMin) || (Obs > ObsMax))
                Obs = NaN;
            end
            VectorData.Lon(cnt) = Lon;
            VectorData.Lat(cnt) = Lat;
            VectorData.Obs(cnt) = Obs;
        end
    end
    scatter(VectorData.Lon,VectorData.Lat,12, VectorData.Obs, 'square', 'filled');
    axis([-180 180 -90 90]);
    title(TitleText, 'FontSize', 10);
    xlabel(LonTag);
    ylabel(LatTag);
    colorbar;
end 


% Cleanup outstanding temporary variables
clear('Lat','Lon','m','n','LatBins','LonBins','Obs','cnt','numAttributes','ObsMax','ObsMin');
clear('LonTag','LatTag','varargin','nVarargs','VectorData');

plot = true;

% ***************************** End **************************************
end

