function toTextFile(data,filename, par)

% USE:
%   Sends data to a tab-delimeted text file with columns representing 
%   each vector of data submitted.
%
% INPUT: 
%       data = n x m matrix where each column will be a column of the
%       text file
%       filename = output file name which can include path (string)
%
% OPTIONAL INPUT:
%           par.makeHeader = 'True' will add a header to file ('True'/'False)
%           par.header = 1 x m cell of header titles to be added to file
%           par.saveFormatIn = 'True' generates a formatIn string at the top
%           of the file for easy reading into matlab
% OUTPUT: 
%       text file with title and location defined by filename. Default path
%       is current directory

%generate file for writing
fileID = fopen(filename,'w');


%format data for printing to file
data = data';
[m,~] = size(data);
%make header - optional
if exist('par','var') && isfield(par, 'makeHeader') && strcmp(par.makeHeader, 'True')
    headerFormatIn = '';
    formatIn = '';
    for  ii = 1:length(par.header)
        headerFormatIn = [headerFormatIn '%s\t'];  
    end
    
    for ii = 1:m
         if isa(data(ii,1),'float')
             formatIn = [formatIn '%f\t'];
         elseif isa(data(ii,1),'integer')
             formatIn = [formatIn '%d\t'];
         else
             formatIn = [formatIn '%s\t'];
         end
    end
    headerFormatIn = [headerFormatIn '\n'];
    par.header{1,1} = ['#' par.header{1,1}];

    if exist('par','var') && strcmp(par.saveFormatIn, 'True')
        fprintf(fileID,'%s\n',['#Data Format: ' formatIn ]);
    end
    fprintf(fileID,headerFormatIn,par.header{1,:});
end

%print to file
fprintf(fileID,[formatIn '\n'],data);

%close file
fclose(fileID);
end