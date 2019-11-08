function pdf = histo(data,plotpar)
% Info:
%   Creates and plots the histogram (pdf) 

% Input:
%   data = row vector containing data to create the histogram with uniform bins 
%   normalization = factor to divide all data by
%       'std' = standard deviation of dataset
%       1 = default
%   bindwidth = integer, distance between bins
%   plotpar = type of marker to use for data points
%       '' = 
%   holdon = 'true' if plotting more than one set on the same graph,
%   default is 'false'

%Output:
%   plot of histogram of input data

%% 

[N, edges] = histcounts(data,'normalization','pdf');

% center bins
for jj = 1:length(edges)-1
    bins(jj) = (edges(jj)+edges(jj+1))/2;     
end

pdf.bins = bins;
pdf.count = N;
pdf.sigma = nanstd(data(:));
pdf.mu = nanmean(data(:));

%plot
figure()
scatter(pdf.bins,pdf.count,plotpar.markersize,plotpar.marker, 'markeredgecolor',plotpar.color)
xlim([min(data(:)),max(data(:))])
t = title('Normalized Probability Distribution');
x = xlabel(plotpar.xlabel);
y = ylabel(plotpar.ylabel);

set(t,'Interpreter','Latex')
set(x,'Interpreter','Latex')
set(y,'Interpreter','Latex')
set(gca,'fontsize',plotpar.fontsize)
set(gcf, 'Position', [45, 1000000, 1200, 900]);           % [left bottom width height]
set(gcf,'color','w');

if plotpar.holdon =='true'
    hold on;
else
    hold off;
end

end