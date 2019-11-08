# plots structure functions of order 1,2,3,4, and 6.

# INPUTS:
#   l
#   t
#   lci
#   tci
#   par: dictionary with key parameters for plotting
#       order
#       


# OUTPUTS:
#   handle of plot figure


def VelocitySF_plot(sf,r,lci,uci,par):
    
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import colors as mcolors
    sys.path.insert(0, '/Users/jennapearson/Google Drive (jenna_pearson@brown.edu)/Functions')
    from figprops import figprops

    # set default figure properties
    figprops()

    fig=plt.figure()
    ax =fig.add_subplot(111)
    cicolor = par['cicolor']
    lncolor = par['lncolor']
    lnstyle = par['lnstyle']
    titleheight = 1.04
    alpha = par['alpha']

    # structure functions
    ax.loglog(r,tsf,color = lncolor, linestyle = lnstyle)
    ax.fill_between(r, uci, lci,facecolor=cicolor, edgecolor = cicolor, alpha=alpha)

    ax.set_title('Second Order Velocity Structure Function', y=titleheight)
    ax.set_xlabel('$r \ (km)$')
    ax.set_ylabel('$(m/s)^2$')
    ax.set_xlim([min(r), max(r)])
    ax.set_ylim([min(sf)-0.001, max(sf) + 0.001])
    

    # guideline plotting
    glfontsz = par['glfontsz']
    glcoef = par['glcoef']
    glcolor = par['glcolor']
    x = np.arange(1,5)
    y1 = np.power(x,(2./3))
    y2 = np.power(x,2)
    y3 = x

    ax.loglog(x,glcoef*y1, color = glcolor)
    ax.loglog(x,glcoef*y2, color =glcolor)
    ax.loglog(x,glcoef*y3, color =glcolor)
    
    # guideline labels
    ax.text(0.33, 0.95,'$r^2$',
        fontsize = glfontsz,
        horizontalalignment='center',
        verticalalignment='center',
        transform = ax.transAxes)
    ax.text(0.335, 0.75,'$r$',
        fontsize = glfontsz,
        horizontalalignment='center',
        verticalalignment='center',
        transform = ax.transAxes)
    ax.text(0.345, 0.69,'$r^{2/3}$',
        fontsize = glfontsz,
        horizontalalignment='center',
        verticalalignment='center',
        transform = ax.transAxes)
                

    # aesthetics
    plt.grid(True,which="both",ls=":")
                 
    plt.show()

    return(fig,ax)
