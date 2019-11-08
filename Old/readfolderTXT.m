function data = readTXTfolder(path,filetype,formatspec)

inputdata = [path '*.' filetype];
datanames = dir(inputdata); % get list of all .nc files in directory
datanames = struct2cell(datanames);
datanames = natsortfiles(datanames(1,:));
nf = length(datanames);
for dd = 1:nf
    input = fopen([path datanames{1}], 'r');          %   define the data types for each column:f=float
    tline = fgetl(input);                                         %   gets next line
    while ~isempty(tline) && any(strncmp(tline, {'*', '#', '/'},1))
       filepos = ftell(input);                                    %   returns current position in the file
       tline = fgetl(input);                                      %   gets next line
    end

    fseek(input, filepos, 'bof'); 
    if nargin<4
        D = textscan(input, formatspec);                            %   read file
    else
        D = textscan(input, formatspec, 'Headerlines',hdrlns);  
    end

    fclose(input); 
end
end