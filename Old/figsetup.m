function figsetup(clr)
    % groot work
    if nargin>=1
        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor', 'DefaultTextcolor'},{clr,clr,clr,clr})
    end
    set(groot,'defaulttextinterpreter','latex');  
    set(groot, 'defaultAxesTickLabelInterpreter','latex');  
    set(groot, 'defaultLegendInterpreter','latex'); 
end