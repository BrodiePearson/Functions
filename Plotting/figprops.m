function figprops(gcf)

    %figure properties
    %set(gcf, 'Position', [45, 1000000, 1300, 1000]);           % for big screen [left bottom width height]
    set(gcf, 'Position', [45, 1000000, 1000, 700]); % for laptop screen
    set(gcf,'color','w');
    
    %axes properties
    set(findall(gcf,'-property','FontSize'),'FontSize',30)
    axesHandles = get(gcf,'children');
    axesHandles = findall(0,'type','axes');
    %set(axesHandles, 'Fontsize', 30);
    set(groot,'defaulttextinterpreter','latex');  
    set(groot, 'defaultAxesTickLabelInterpreter','latex');  
    set(groot, 'defaultLegendInterpreter','latex');   
    
    %line properties
    lineHandles = findall(0,'type','line');
    set(lineHandles, 'linewidth', 3); 
end
    