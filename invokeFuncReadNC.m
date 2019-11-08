 % ***********  Call Script and Example Usage of generic netCDF Read Functions  ***********
%  
%  SCRIPT: invokeFuncReadNC.m
%  
%  DESCRIPTION
%   - Illustrates how readNC functions are called from within a simple Matlab script
%   - This script can be used to read generic netCDF data files, display metadata & data listings.
%
%  USAGE
%   The script automatically invokes the readNC.m function to read a
%   generic netCDF datafile and return complete metadata and variable data.
%
%  DEPENDENCIES
%   The script invokes the following CAP reader functions and source files written in Matlab and provided by PO.DAAC:
%       readNC.m   readNC()
%   Note: see the source .m files for complete function documentation.
%
%  INPUTS   ("User Inputs" code block section)
%       fpath: file path
%       fname: file name of source  NetCDF data file 
%       DisplayListings: True/False setting to output data listings sequentially by Dataset parameter
%            
%  OUTPUTS
%       globalDimensions: structure containing global dimension attributes for arrays
%       globalMetadata: structure containing global file attributes
%       VarData: multidimensional data structure (dimensioned on #variables)
%       optional Listings of dataset attributes & data to screen
%  
%  NOTES
%      1. This read software was created using Matlab version 7.14
%      2. Please email all comments and questions concerning these routines
%           to podaac@podaac.jpl.nasa.gov.
%
%  CREATED:
%       4/5/2016: Vardis Tsontos, PO.DAAC, NASA-JPL CalTech
%
%  UPDATED:
%       .
%     
%
%======================================================================
% Copyright (c) 2016, California Institute of Technology
%======================================================================



% ****** User Inputs Section ***********            % User inputs entered here
Path = 'C:\Users\vtsontos\Desktop\';                % Enter the data directory path here
Name = 'SSS_OI_7D_20112392011245_V40.nc';           % Enter the data filename to be read
ListFlag = false;                                   % Output data listings: specify true or false

% ---------------------------------------------------------------------------------------------------------------------
  


% *******************************    Main   *******************************

% executes the function and passes defined attribute values)

    [globalDimensions, globalMetadata, VarData]=readNC(Path,Name,ListFlag);

%--------------------------------------------------------------------------
