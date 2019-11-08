function [numGlobalAttributes,numDatasets,numGroups,globalMetadata,dataset,group,groupInfo] = readHDF5(fPath, fName, ListFlag)
%  
%  FUNCTION: readHDF5(fpath,fname,ListFlag)
%
%  DESCRIPTION
%   This function reads into memory in its entirity any arbitrary HDF5 data file selected by the user, with an option for
%   additionally displaying file contents as listings to screen.  This includes all available hierarchically structured 
%   file metadata and data elements, including global attributes and both Dataset and Group level metadata.  Input 
%   arguments require specification of the file path/name and listing output flag.
%
%  USAGE
%   Call the "readHDF5(fpath,fname,ListFlag)" function from either the Matlab command
%   line or within a script using suitable arguement values
%
%  INPUTS ARGUMENTS      readHDF5(fPath,fName,ListFlag)
%       fpath: file path (inlude terminal "\" in the path string   eg. C:\data\ )
%       fname: name of source Aquarius HDF5 L2 or L3 data file
%       ListFlag: True/False flag setting to output data listings sequentially by hierarchical structure and element
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



FileName = [fPath fName];

%Read and store Global level HDF5 File Metadata
FileInfo = h5info(FileName);
numGlobalAttributes = length(FileInfo.Attributes);
if (numGlobalAttributes > 0 )
    globalAttributes = FileInfo.Attributes;
    for n =1: numGlobalAttributes
        if n == 43
            n=n+1;
        end
        if isnumeric(globalAttributes(n).Value)
            attribVal = num2str(globalAttributes(n).Value(:)');
        else
            attribVal = globalAttributes(n).Value;
        end    
        globalMetadata(n).Attribute = [globalAttributes(n).Name  ':  ' attribVal];
    end
    clear('n','globalAttributes','attribVal');    
end


%Read and store available Dataset level HDF5 File Metadata & Data
numDatasets = length(FileInfo.Datasets);    
if (numDatasets > 0 )
    datasetAttributes = FileInfo.Datasets;
    for n =1: numDatasets
        numDatasetAttributes = length(FileInfo.Datasets(n).Attributes);
        Temp = '';
        for m =1: numDatasetAttributes
            if isnumeric(datasetAttributes(n,1).Attributes(m).Value)
                attribVal = num2str(datasetAttributes(n,1).Attributes(m).Value(:)');
            else
                attribVal = datasetAttributes(n,1).Attributes(m,1).Value;
            end            
            Temp = [Temp ' , ' datasetAttributes(n,1).Attributes(m,1).Name  ':  ' attribVal];
        end
    dataset(n).Name = datasetAttributes(n).Name;
    dataset(n).numAttributes = numDatasetAttributes;
    dataset(n).Attributes = strtrim(Temp);
    dataset(n).Data = h5read(FileName,['/' datasetAttributes(n).Name]);    
    end
    clear('n','m','numDatasetAttributes','datasetAttributes','attribVal','Temp');
else
    dataset =0;
end


%Read and store available Group level HDF5 File Metadata & Data
groupInfo = FileInfo.Groups; 
numGroups = length(groupInfo);
if (numGroups > 0 )
    for n =1: numGroups
      %Get Group Level attributes
        numGroupAttributes = length(groupInfo(n,1).Attributes);
        Temp = '';
        if (numGroupAttributes > 0 )
            for m =1: numGroupAttributes
                if isnumeric(groupInfo(n,1).Attributes(m).Value)
                    attribVal = num2str(groupInfo(n,1).Attributes(m).Value(:)');
                else
                    attribVal = groupInfo(n,1).Attributes(m,1).Value;
                end            
                Temp = [Temp ' , ' groupInfo(n,1).Attributes(m,1).Name  ':  ' attribVal];
            end
        end
        numGroupDatasets = length(groupInfo(n,1).Datasets);
        group(n).Name = groupInfo(n).Name;
        group(n).numAttributes = numGroupAttributes;
        group(n).Attributes = strtrim(Temp);        
        group(n).numDatasets = numGroupDatasets;        

     %Get Group Dataset Level Metadata and Data
        for i =1: numGroupDatasets
            group(n).Dataset(i).Name = groupInfo(n,1).Datasets(i,1).Name;
            group(n).Dataset(i).Data = h5read(FileName, [group(n).Name '/' group(n).Dataset(i).Name]);
            numGroupDatasetAttributes = length(groupInfo(n,1).Datasets(i,1).Attributes);
            group(n).Dataset(i).numAttributes = numGroupDatasetAttributes;
            if (numGroupDatasetAttributes > 0 )
                Temp = '';
                for j =1: numGroupDatasetAttributes
                    if isnumeric(groupInfo(n,1).Datasets(i,1).Attributes(j,1).Value)
                        attribVal = num2str(groupInfo(n,1).Datasets(i,1).Attributes(j,1).Value(:)');
                    else
                        attribVal = groupInfo(n,1).Datasets(i,1).Attributes(j,1).Value;
                    end            
                    Temp = [Temp ' , ' groupInfo(n,1).Datasets(i,1).Attributes(j,1).Name  ':  ' attribVal];
                end
                group(n).Dataset(i).Attributes = strtrim(Temp);                
            end
        end
    end
    clear('i','j','n','m','numGroupAttributes','numGroupDatasets','numGroupDatasetAttributes','attribVal','Temp');
else
    group =0;
end

clear('FileInfo', 'fPath','fName');



% --- Displays listings of all HDF5 file global attributes, Datasets & Group data and metadata 
%           if the user has selected to do so (DisplayListings flag = True) ---

if ListFlag
    disp(' ');
    disp(['File Name: ' FileName]);
    disp(['Number of Global Attributes: ' num2str(numGlobalAttributes)]);
    disp(['Number of Datasets: ' num2str(numDatasets)]);
    disp(['Number of Groups: ' num2str(numGroups)]);
    disp(' ');    
    disp('     Global Attributes');
    disp('------------------------------------'); 
    if (numGlobalAttributes > 0)
        for n =1:numGlobalAttributes
            disp(num2str(globalMetadata(n).Attribute));
        end
    end;    
    if (numDatasets > 0)
        disp(' ');
        disp('--------- Press Space to Continue with Dataset Attributes & associated Data listings  --------------');
        pause;
        disp(' ');
        disp('     Dataset Attributes & Data');
        disp('------------------------------------');
        for n = 1:numDatasets
            disp(['Dataset # ' num2str(n)]);
            disp(['Dataset Name: ' dataset(n).Name]);
            disp(['Dataset Attributes: ' dataset(n).Attributes]);
            disp('Dataset Values: ');
            disp(dataset(n).Data);
            if n < numDatasets
                disp(['--------- Press Space to Continue for Remaining ' num2str(numDatasets-n) ' Datasets  --------------']);
                pause;
            end       
        end
    end
    if (numGroups > 0)
        disp(' ');
        disp('--------- Press Space to Continue with Group Attributes & associated Data listings  --------------');
        pause;
        disp(' ');
        disp('     Group Attributes & Data');
        disp('------------------------------------');
        for n = 1:numGroups
            disp(['Group # ' num2str(n)]);
            disp([num2str(n) ' Group Name: ' group(n).Name]);
            disp([num2str(n) ' Group Attributes: ' group(n).Attributes]);
            disp([num2str(n) ' # Group Datasets: ' num2str(group(n).numDatasets)]);
            
            for m = 1:group(n).numDatasets
                disp([num2str(n) ' Dataset #: ' num2str(m)]);
                disp([num2str(n) ':' num2str(m) ' Dataset Name: ' group(n).Dataset(m).Name]);
                disp([num2str(n) ':' num2str(m) ' Dataset Attributes: ' group(n).Dataset(m).Attributes]);
                disp([num2str(n) ':' num2str(m) ' Dataset Data: ' ]);
                disp(group(n).Dataset(m).Data);
                
                if m < group(n).numDatasets
                    disp(['--------- Press Space to Continue for Remaining ' num2str(group(n).numDatasets-m) ' Datasets of Group ' num2str(n)  ' --------------']);
                    pause;
                end
            end
        end    
    end
end

% ***************************** End **************************************
