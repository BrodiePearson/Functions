function superxlabel(gcf,title,fontsz,color)
 
if nargin < 3
    fontsz = 20;
    color = 'k';
elseif nargin < 4
    color = 'k';
end
    
%axes properties
    axesHandles = get(gcf,'children');
    axesHandles = findall(0,'type','axes');
    pos = get(axesHandles, 'Position');
    for ii = 1:size(pos,1)
        pos{ii}(1,2) = 0.11;
        set(axesHandles(ii), 'Position', pos{ii})
    end
   
annotation('textarrow', [0.5 0.5],[0.07 0.5], ...
    'HeadStyle','none', ...
    'LineStyle', 'none', ...
    'String', title, ...
    'Fontsize', fontsz, ...
    'Color', color, ...
    'HorizontalAlignment', 'center',...
    'Interpreter', 'latex')
end