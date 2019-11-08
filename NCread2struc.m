function [VarData,globalDimensions, globalMetadata] = NCread2struc(fPath, fName, flipdim)

%  INFO:
%   This function reads an arbitrary user defined NetCDF file and returns a data structure array (VarData) 
%   containing complete variable metadata and data. It also returns
%   variables containing information on source file global metadata (globalMetadata) and variable dimension (globalDimensions).
%   Input arguements require file path and name specification, and provide the user with options to dump both
%   dataset attributes & data as listings to screen for review.
%
%  INPUTS:
%       fname: file name of source NetCDF data file
%       fpath: file path
%       DisplayListings: True/False flag setting to output data listings sequentially by Dataset parameter
%
%  OUTPUTS:
%       globalMetadata: structure containing global file attributes
%       VarData: multidimensional data structure (dimensioned on number of variables)
%           Name: parameter/variable name
%           Attributes: concatenated parameter attributes (types & values)
%       optional Listings of dataset attributes & data to screen



% Read file Dimensions metadata into GlobalDims structure
FileName = [fPath fName];
FileInfo = ncinfo(FileName);
FileFormat = FileInfo.Format;
globalDimensions = FileInfo.Dimensions;
numGlobalDimensions = length(FileInfo.Dimensions);
for n=1: numGlobalDimensions 
    globalDimensions(n).Attribute = [globalDimensions(n).Name  ': ' num2str(globalDimensions(n).Length)];
end
    
% Read file metadata & global netCDF file attributes into globalMetadata structure
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


% Read NetCDF file data and variable attributes into VarData structure
for n = 1:numDatasets
    
    if flipdim
        numAttributes = length(FileInfo.Variables(n).Attributes);
        Temp ='';
        for m = numAttributes:-1:1
            Temp = [Temp '   ' FileInfo.Variables(n).Attributes(m).Name  ': ' num2str(FileInfo.Variables(n).Attributes(m).Value)];
        end
        VarData(n).Name = FileInfo.Variables(n).Name;
        VarData(n).Attributes = strtrim(Temp);
        tempvardata = ncread(FileName,['/' FileInfo.Variables(n).Name]);
        VarData(n).Data = tempvardata';
    else
        numAttributes = length(FileInfo.Variables(n).Attributes);
        Temp ='';
        for m = numAttributes:-1:1
            Temp = [Temp '   ' FileInfo.Variables(n).Attributes(m).Name  ': ' num2str(FileInfo.Variables(n).Attributes(m).Value)];
        end
        VarData(n).Name = FileInfo.Variables(n).Name;
        VarData(n).Attributes = strtrim(Temp);
        VarData(n).Data = ncread(FileName,['/' FileInfo.Variables(n).Name]);
    end
end


% Clear Temp variables
clear('FileInfo', 'fPath', 'fName','globalAttributes');
clear('m', 'n','Temp');
clear('numAttributes','attribVal');

