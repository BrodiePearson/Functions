function cellboxplot(data)
% data is a cell containing one array for each boxplot to be added.
% go through lengths of cells to define groups 
grps = [];
grpdata = [];
for cc = 1:length(data)
    grpdata = [grpdata data{cc}'];
    grps = [grps cc.*ones(1,length(data{cc}))];
end

hold on;
boxplot(grpdata,grps,'boxstyle','filled')

end