import pickle
import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import UnivariateSpline

num_tot=5

cm = plt.get_cmap('viridis')
colors = [cm(1.4*float(n)/num_tot) for n in np.arange(num_tot)]

for X in ['CDM','SI3','SI10','vdXsec']:
    f,ax=plt.subplots(1,2,figsize=(6,3))
    plt.subplots_adjust(wspace=0)
    f.suptitle(X,fontsize=15)
    for i in [0,1]:
        ax[i].set_ylim([0,1])
        ax[i].set_xlim([0,1])
        ax[i].set_yticks(np.linspace(0,1,6))
    ax[0].set_xticks(np.linspace(0,.8,5))
    ax[1].set_xticks(np.linspace(0,1,6))
    ax[0].set_ylabel('b/a',fontsize=15)
    ax[0].set_xlabel(r'R/R$_{vir}$',fontsize=15)
    ax[1].set_ylabel('c/a',fontsize=15)
    ax[1].yaxis.set_label_position("right")
    ax[1].yaxis.tick_right()
    ax[1].set_xlabel('R [kpc]',fontsize=15)

    data = pickle.load(open(f'../DataFiles/DarkMatterShapes.Top200.{X}.pickle','rb'))
    num_plotted,halo = 0,1
    while num_plotted < num_tot:
        try:
            h = data[str(halo)]
            r = h['rbins']/h['rbins'][-1]
            sm_b = UnivariateSpline(r,h['b_pyn'],k=3)
            sm_c = UnivariateSpline(r,h['c_pyn'],k=3)
            ax[0].plot(r,sm_b(r),color=colors[num_plotted])
            ax[1].plot(r,sm_c(r),color=colors[num_plotted])
            num_plotted+=1
            halo+=1
        except:
            halo+=1
    f.savefig(f'../Plots/ShapeProfiles.{X}.png',bbox_inches='tight',pad_inches=.1)