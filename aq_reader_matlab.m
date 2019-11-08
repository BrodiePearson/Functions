function [lons,lats,sss] = aq_reader_matlab(fn)
%
%  FILENAME: aq_H5_reader_podaac.m
%
%  BATCH USAGE: matlab < /path/to/this/file/aq_H5_reader_podaac(fn)
%  
%  USAGE: aq_H5_reader_podaac(fn)
%
%  User options can be changed below by changing flags from true/false
%    to false/true. These options involve displaying or saving metadata,
%    and plotting options.      
%
%  INPUT
%  fn: filename of hdf file w/path of aquarius hdf5 file
%
%  OUTPUT
%     Level 3
%        lons: latitudes            360x8byte
%        lats: longitudes           180x8byte
%        sss:  sea surface salinity 360x180x4byte 
%     Level 2
%        lons: latitudes            3x4084x4byte %middle dimmension will vary
%        lats: longitudes           3x4084x4byte
%        sss:  sea surface salinity 3x4084x4byte
%     optional figure to screen
%     optional attributes file
%     optional attributes to screen
%  
%  DESCRIPTION
%  Level 3
%     This file reads Level 3 hdf5 aquarius products and returns a 2 dimensional
%     Sea Surface salinity variable, and longitude and latitude vectors.
%  Level 2
%     This file reads Level 2 hdf5 aquarius products and returns 2 dimensional 
%     latitude, longitude, and sea surface salinity variables.  One dimension is
%     3, 1 for each radiometer and the other for the number of observations.
%  Levels 2 & 3 all can optionally display, or print to file the attributes of 
%     the granule as well as plot the data.
%
%  NOTES:
%
%      1. This read software was created using Matlab version 7.4.0.287.
%
%      2. Please email all comments and questions concerning these routines
%           to podaac@podaac.jpl.nasa.gov.
%
%  CREATED:
%       2010 October 21: G. Foti
%
%======================================================================
% Copyright (c) 2010, California Institute of Technology
%======================================================================
%
[pathstr,fsname,ext] = fileparts(fn);
fsname = [fsname ext]; % file name w/o direcory
if findstr(fsname,'L3')
   lev3 = true;
else
   lev3 = false;
end

if lev3 & findstr(fn,'main')
    error('This reader only reads the along-track L2 data and the gridded L3 data. It does not read the L3 binned data');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER OPTIONS: 

%using diary for this.  could use fprintf but would add lines to code 
attr_file = false;%to write attributes to file set to true

attr_disp = true; %write attributes to screen or set to false
make_fig  = true; %plot data or set to false
examine_params = false; %pauses function to examine all parameters

%for level 3 ignore sss values outside fig_filt range when plotting
%for level 2 set range of colorbar to these fig_filt values when plotting SSS
fig_filt = [28 40]; %lower salinity must be first

%BELOW only for plotting level 2 data
if make_fig & ~lev3    %if plotting L2 file which parameters to plot?
   pparams = {'SSS'}  %only plot SSS parameter
   %pparams = {'all'}  %plot all parameters

   %example of how to customize which of the 79 L2 parameters to plot
   %pparams = {'SSS' 'scat_wind_speed' 'anc_surface_temp'}

   display_time = 0.5; %how many seconds should each parameter be displayed
                       %set to 0 to keep plot displayed until return is entered
end
%END USER OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIPLAY ACKNOWLEDGEMENT INFO
%
disp(' ')
disp('If this data product is used for publication please include the')
disp('following acknowledgement:')
disp(' ')
disp('''The Aquarius Sea Surface Salinity data was obtained from the JPL')
disp('Physical Oceanography Distributed Active Archive Center (PODAAC).''')
disp(' ')
disp('and if appropriate:')
disp(' ')
disp('''Sea surface salinity image courtesy of Physical Oceanography')
disp('Distributed Active Archive Center (PODAAC).''')
disp(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if capturing metadata info in a file use diary function.
if attr_file
   attr_disp = true;
   diary_nam = [fsname '_att'];
   if exist(diary_nam,'file')
      system(['rm ' diary_nam]);
   end
   diary([fsname '_att'])
end
 
if lev3
   %lats and lons on 1 degree boundary
   lats = -89.5:89.5;
   lons = -179.5:179.5;
else
   beam_clat = hdf5read(fn,'/Navigation/beam_clat');
   beam_clon = hdf5read(fn,'/Navigation/beam_clon');
end

%Open HDF file
fileinfo = hdf5info(fn);
toplevel = fileinfo.GroupHierarchy;
numAttributes = length(toplevel.Attributes);
if attr_disp, disp(['number of Global Attributes = ' num2str(numAttributes)]), end 

if attr_disp % capturing attribute info?
   for i = 1:numAttributes
      if attr_disp
         disp(['Attribute ' num2str(i)])
         disp(['Attribute Name-> '  toplevel.Attributes(i).Name])
         if        ischar(toplevel.Attributes(i).Value)
            disp(['Attribute Value->' toplevel.Attributes(i).Value])
         elseif isnumeric(toplevel.Attributes(i).Value)
            disp(['Attribute Value->' num2str(toplevel.Attributes(i).Value(:)')])
         else
            disp(['Attribute Value->' toplevel.Attributes(i).Value.data])
         end
      end
   end
end

numGroups = length(toplevel.Groups);
if attr_disp
    display(' ')
    display(['File contains ' num2str(numGroups) ' groups'])
end

for i = 1:numGroups
   numGroupDsets = length(toplevel.Groups(i).Datasets);
   if attr_disp
      display(' ')
      disp(['Group ' num2str(i) ' contains ' num2str(numGroupDsets) ' datasets'])
      disp(['Group Name-> '  toplevel.Groups(i).Name])
   end

   for j = 1:numGroupDsets
      if attr_disp, disp(['GrpDset ' num2str(j)]), end
      [dum,grpdsname,dum] = fileparts(toplevel.Groups(i).Datasets(j).Name);
      if attr_disp, disp(['GrpDset Name->' grpdsname]), end
      eval([grpdsname ' = hdf5read(fn, toplevel.Groups(i).Datasets(j).Name);']);
      if attr_disp
          disp(['GrpDataset Dims-> ' ...
                num2str(toplevel.Groups(i).Datasets(j).Dims)])
      end

      numGrpDatasetAttributes=length(toplevel.Groups(i).Datasets(j).Attributes);
      if attr_disp
         disp(['number of dataset attributes-> ' ...
                num2str(numGrpDatasetAttributes)])
      end

      for k = 1:numGrpDatasetAttributes
         if attr_disp, disp(['   GrpDset Attribute ' num2str(k)]), end
         [dum,grpdsatname,dum]= ...
            fileparts(toplevel.Groups(i).Datasets(j).Attributes(k).Name);
         if attr_disp, disp(['   GrpDset Att. Name->' grpdsatname]), end

         % get slope, intercept and Units
         if strcmp(grpdsatname,'Slope')
            slope = toplevel.Groups(i).Datasets(j).Attributes(k).Value;
         elseif strcmp(grpdsatname,'Intercept')
            intercept = toplevel.Groups(i).Datasets(j).Attributes(k).Value;
         elseif strcmp(grpdsatname,'units')
            units = toplevel.Groups(i).Datasets(j).Attributes(k).Value.Data;
         elseif strcmp(grpdsatname,'long_name')
            longname = toplevel.Groups(i).Datasets(j).Attributes(k).Value.Data;
         end

         if attr_disp 
            try
            disp(['   GrpDset Att. Value->' num2str(toplevel.Groups(i).Datasets(j)...
                                     .Attributes(k).Value)])
            catch
               %must have more than 1 element
               try
                  disp(['   GrpDset Att. ValueName->' toplevel.Groups(i).Datasets(j)...
                                           .Attributes(k).Value.Name])
                  disp(['   GrpDset Att. ValueDATA->' num2str(toplevel.Groups(i)...
                                         .Datasets(j).Attributes(k).Value.Data)])
               catch
                  disp(['   GrpDset Att. Value->' num2str(toplevel.Groups(i).Datasets(j)...
                                           .Attributes(k).Value')])
               end
            end
         end
      end
      if make_fig & ~lev3
         plot_this = false;
         for pparam = pparams
            if strcmp(pparam,'all') | strcmp(pparam,grpdsname)
               plot_this = true;
               break
            end
         end

         if plot_this
            dssize = eval(['size(' grpdsname ');']);
            if length(dssize) == 2 & dssize(1) == 3
                  titout = [longname ' (' units ')'];
                  eval(['parplot = ' grpdsname ';']);

                  if strcmp(grpdsname,'SSS') | strcmp(grpdsname,'scat_wind_speed')...
                   | strcmp(grpdsname,'anc_surface_temp')
                     bad_flag = find(parplot < 0);
                     parplot(bad_flag) = NaN;
                  end

                  scatter(beam_clon(1,:),beam_clat(1,:),5,parplot(1,:))
                  hold on
                  scatter(beam_clon(2,:),beam_clat(2,:),5,parplot(2,:))
                  scatter(beam_clon(3,:),beam_clat(3,:),5,parplot(3,:))
                  if strcmp(grpdsname,'SSS')
                     caxis([fig_filt(1) fig_filt(2)])
                  end
                  colorbar
                  grid on
                  hold off
                  grpdsname = strrep(grpdsname,'_',' '); 
                  title(titout)
                  if ~display_time
                     display('Enter ''return'' to continue')
                     keyboard
                  else
                     pause(display_time)
                  end
            end 
         end 
      end
   end
end

numDatasets = length(toplevel.Datasets);
if attr_disp, disp(['num datasets = ' num2str(numDatasets)]), end
for i = 1:numDatasets
   [dum,dsname,dum] = fileparts(toplevel.Datasets(i).Name);
   if attr_disp
      disp(['Dataset ' num2str(i)])
      disp(['Dataset Name-> '  dsname])
   end
   eval([dsname ' = hdf5read(fn, dsname);']);
   if attr_disp
      disp(['Dataset Dims-> '  num2str(toplevel.Datasets(i).Dims)])
   end

   numDatasetAttributes = length(toplevel.Datasets(i).Attributes);
   if attr_disp
      disp(['num Dset Attributes = ' num2str(numDatasetAttributes)])
   end

   for j = 1:numDatasetAttributes
      [dum,dsatname,dum]=fileparts(toplevel.Datasets(i).Attributes(j).Name);
      if attr_disp
         disp(['   Dset Attribute' num2str(j)])
         disp(['   Dset Att. Name->' dsatname])
      end 

      if strcmp(dsatname,'Slope')
         slope = toplevel.Datasets(i).Attributes(j).Value;
      elseif strcmp(dsatname,'Intercept')
         intercept = toplevel.Datasets(i).Attributes(j).Value;
      end

      if attr_disp 
         try
            disp(['   Dset Att. Value->' num2str(toplevel.Datasets(i)...
                                     .Attributes(j).Value)])
         catch
            %must have more than 1 element
            disp('   dataset attribute contains more than 1 element')
            disp(['   Dset Att. ValueName->' toplevel.Datasets(i)...
                                  .Attributes(j).Value.Name])
            disp(['   Dset Att. ValueDATA->' num2str(toplevel.Datasets(i)...
                                  .Attributes(j).Value.Data)])
         end
      end
   end
   if strcmp(dsname,'l3m_data') & lev3
      sss = single(l3m_data) * slope + intercept;
      [rows cols] = size(sss);
      if make_fig

         if strfind(fn,'wind')
            noplot = find(sss < 0);
         else
            noplot = find(sss <= fig_filt(1) | sss >= fig_filt(2) | sss < 0);
         end

         plotdat = double(sss); 
         plotdat(noplot) = NaN;
  
         pcolor(lons,lats,plotdat(:,cols:-1:1)')
         shading flat
         grid on
         colorbar
         fsname = strrep(fsname,'_','\_'); 
         title(fsname)
      end
   end
end

if lev3 & ~exist('sss', 'var')
    allsss=h5read(fn,[toplevel.Groups(1).Name '/SSS']);
    sss=allsss.SSS_sum;
    [rows cols] = size(sss);
    if make_fig
        
        if strfind(fn,'wind')
            noplot = find(sss < 0);
        else
            noplot = find(sss <= fig_filt(1) | sss >= fig_filt(2) | sss < 0);
        end
        
        plotdat = double(sss);
        plotdat(noplot) = NaN;
        
        pcolor(lons,lats,plotdat(:,cols:-1:1)')
        shading flat
        grid on
        colorbar
        fsname = strrep(fsname,'_','\_');
        title(fsname)
    end
end
   
if ~lev3
   sss = SSS;
   lons = beam_clon;
   lats = beam_clat;
end

if attr_file 
   diary off
end
if examine_params
   display('Enter ''return'' to exit function')
   keyboard
end
