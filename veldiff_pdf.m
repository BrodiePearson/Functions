clear all; close all; clc;

%% set paths and load data
% 
% % LASER Xband Radar
% launch = 'Xband_Radar';
% path = '/Users/jennapearson/Google Drive (jenna_pearson@brown.edu)/Baylor & Jenna/Projects/Eulerian-Lagrangian SFs (Obs)/LASER/Data/';
% dataname = 'XBandRadar_LASER_processed_and_formatted.mat';
% 
% load([path 'Processed and Formatted/' dataname])


%LASER DRIFTERS
launch = 'All';
path = '/Users/jennapearson/Google Drive (jenna_pearson@brown.edu)/Baylor & Jenna/Projects/Eulerian-Lagrangian SFs (Obs)/LASER/Data/';
dataname = ['LASER_' launch '_SPOT_Drifters_processed_and_formatted.mat'];

load([path 'Processed and Formatted/' dataname])
 
% % load GLAD S1 data as a test
% path = '/Users/jennapearson/Google Drive (jenna_pearson@brown.edu)/Research /Gulf of Mexico/2012/Drifters/GLAD/CODE/Data/Processed/';
% dataname = 'GLAD_S1.mat';
% 
% load([path dataname])
%%  Get relative velocities for each time step
[ul,ut,r] = rel_vel_short(lon,lat,u,v);

%combine all time steps
ultemp = ul(:);
uttemp = ut(:);
rtemp = r(:)/1000;

%% Bin by scale

% Generate Bin Vector
rBins = 0:.25:ceil(max(rtemp));

% Bin Data
for ii = 1:length(rBins)
    disp(ii)
    %find indicies for bins
    ind = find(abs(rtemp - rBins(ii)) < 0.1);  
    ulbinned{ii} = ultemp(ind);                           
    utbinned{ii} = uttemp(ind);                         
end


% Save data

struc.ul = ul;
struc.ut = ut;
struc.r = r/1000;
struc.rBins = rBins;
struc.ulbinned = ulbinned;
struc.utbinned = utbinned;

outdataname = ['LASER_' launch '_SPOT_Drifters_PDF.mat'];
%outdataname = 'LASER_Xband_Radar_PDF.mat';

save([path 'PDFs/' outdataname],'-struct','struc','-v7.3');


%%Make a histogram at a given scale
close all;

% small scales
subplot(1,2,1)
% 0.25 km
data = ulbinned{2}/std(ulbinned{2});
empdf = fitdist(data,'Kernel');
x = min(data):.1:max(data);
y = pdf(empdf,x);
plot(x,y, 'k.','MarkerSize',12)
hold on 

% 0.5km
data = ulbinned{3}/std(ulbinned{3});
empdf = fitdist(data,'Kernel');
x = min(data):.1:max(data);
y = pdf(empdf,x);
plot(x,y, 'gs')

% 1km
data = ulbinned{5}/std(ulbinned{5});
empdf = fitdist(data,'Kernel');
x = min(data):.1:max(data);
y = pdf(empdf,x);
plot(x,y, 'r^')

% Gaussian 
data = normrnd(0,0.75,[100000,1]);
empdf = fitdist(data,'Normal');
x = min(data):.1:max(data);
y = pdf(empdf,x);
plot(x,y, '--')

title('Small Scales')
ylabel('PDF')
ylim([.0005,1])
xlim([-4.5,4.5])
set(gca, 'YScale', 'log')
legend('0.25km','0.5km','1km')


% large scales
subplot(1,2,2)
% 5 km
data = ulbinned{21}/std(ulbinned{21});
empdf = fitdist(data,'Kernel');
x = min(data):.1:max(data);
y = pdf(empdf,x);
plot(x,y, 'k.','MarkerSize',12)
hold on 

% 15km
data = ulbinned{61}/std(ulbinned{61});
empdf = fitdist(data,'Kernel');
x = min(data):.1:max(data);
y = pdf(empdf,x);
plot(x,y, 'gs')

% 30km
data = ulbinned{81}/std(ulbinned{81});
empdf = fitdist(data,'Kernel');
x = min(data):.1:max(data);
y = pdf(empdf,x);
plot(x,y, 'r^')

% Gaussian 
data = normrnd(0,.75,[100000,1]);
empdf = fitdist(data,'Normal');
x = min(data):.1:max(data);
y = pdf(empdf,x);
plot(x,y, '--')


title('Large Scales')
superxlabel('$\frac{\Delta u}{\sigma}$')
ylim([.0005,1])
xlim([-4.5,4.5])
legend('5km','7km','10km')
set(gca, 'YScale', 'log')

supertitle(['LASER ' launch])

%% Other Old Code
% Glad Data
%path = '/Users/jennapearson/Google Drive (jenna_pearson@brown.edu)/Research /Gulf of Mexico/2012/Drifters/GLAD/CODE/Data/Processed/';
%dataname = 'GLAD_S1.mat';

% hold on
% % 0.25 km
% data = ulbinned{2}/std(ulbinned{2});
% h1 = histogram(data,'Normalization','pdf');
% for jj = 1:length(h1.BinEdges)-1
%     x1(jj) = (h1.BinEdges(jj)+h1.BinEdges(jj+1))/2;     % Declare the average distance between bins as the representative displacement for a given bin ( centering )
% end
% 
% 
% % 0.5 km
% data = ulbinned{3}/std(ulbinned{3});
% h2 = histogram(data,'Normalization','pdf');
% for jj = 1:length(h2.BinEdges)-1
%     x2(jj) = (h2.BinEdges(jj)+h2.BinEdges(jj+1))/2;     % Declare the average distance between bins as the representative displacement for a given bin ( centering )
% end
% 
% 
% % 1 km
% data = ulbinned{5}/std(ulbinned{5});
% h3 = histogram(data,'Normalization','pdf');
% for jj = 1:length(h3.BinEdges)-1
%     x3(jj) = (h3.BinEdges(jj)+h3.BinEdges(jj+1))/2;     % Declare the average distance between bins as the representative displacement for a given bin ( centering )
% end
% 
% figure()
% subplot(1,2,1)
% plot(x1,h1.Values, 'k.','MarkerSize',12)
% hold on
% plot(x2,h2.Values, 'ks')
% plot(x3,h3.Values, 'k^')
% % Gaussian 
% data = normrnd(0,0.5,[100000,1]);
% empdf = fitdist(data,'Normal');
% x = min(data):.1:max(data);
% y = pdf(empdf,x);
% plot(x,y, '--')
% 
% title('Probability Density Function Estimate')
% legend('0.25km','0.5km','1km')
% set(gca, 'YScale', 'log')
% ylim([.0005,1])
% xlim([-4.5,4.5])
