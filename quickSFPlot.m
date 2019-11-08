function quickSFPlot(r,l,t,par)

if par.order == 1
    plot(r,l)
    hold on
    plot(r,t)
    legend('D_L','D_T')
    title('First Order Structure Function')
    xlabel('$r (km)$')
    ylabel('$(m/s)$')
    if isfield(par,'CI')
        grey = [0.5 0.5 0.5];
        plot_ci(r, par.lCI,'PatchColor', grey, 'PatchAlpha', 0.35,'LineColor',grey,...
     'Xscale','log','Yscale', 'log');
        removeCIbound(gca)
    end
elseif par.order == 2
    loglog(r,l)
    hold on
    loglog(r,t)
    legend('D_L','D_T')
    title('Second Order Structure Function')
    xlabel('$r (km)$')
    ylabel('$(m/s)^2$')
    
    if isfield(par,'CI')
        grey = [0.5 0.5 0.5];
        plot_ci(r, par.lCI,'PatchColor', grey, 'PatchAlpha', 0.35,'LineColor',grey,...
     'Xscale','log','Yscale', 'log');
        removeCIbound(gca)
    end
    
    
elseif par.order == 3
    loglog(r,abs(l(l<0)),'o')
    loglog(r,l(l>0),'+')
    hold on
    loglog(r,abs(t(t<0)),'o')
    loglog(r,t(t>0),'+')
    legend('neg D_L','pos D_L','neg D_T','pos D_T')
    title('Third Order Structure Function')
    xlabel('$r (km)$')
    ylabel('$(m/s)^3$')
    
    if isfield(par,'CI')
        grey = [0.5 0.5 0.5];
        plot_ci(r, par.lCI,'PatchColor', grey, 'PatchAlpha', 0.35,'LineColor',grey,...
     'Xscale','log','Yscale', 'log');
        removeCIbound(gca)
    end
end

if isfield(par,'color')
    %line properties
    lineHandles = findall(0,'type','line');
    for ll = 0:length(lineHandles)-1
        set(lineHandles(end-ll), 'color',par.color{ll+1});
    end
end

if isfield(par,'linestyle')
    %line properties
    lineHandles = findall(0,'type','line');
    for ll = 0:length(lineHandles)-1
        set(lineHandles(end-ll), 'linestyle',par.linestyle{ll+1});
    end
    
end

figprops(gcf)
end