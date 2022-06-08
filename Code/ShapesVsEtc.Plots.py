import pickle
import numpy as np
import matplotlib.pylab as plt
from osxmetadata import OSXMetaData

f,ax = plt.subplots(2,2,figsize=(8,8))
plt.subplots_adjust(hspace=0,wspace=0)
for i in [0,1]:
    ax[i][0].set_xlim([7,11])
    ax[i][1].set_xlim([4,100])
    ax[i][1].semilogx()
    for j in [0,1]:
        ax[i][j].set_ylim([0,1])
        ax[i][j].set_xticks([])
        ax[i][j].set_yticks([])
        ax[i][j].tick_params(length=5,labelsize=15)
ax[1][0].set_yticks([0,.2,.4,.6,.8])
ax[0][0].set_yticks([0,.2,.4,.6,.8,1])
ax[1][0].set_xticks([7,8,9,10,11])
ax[1][1].set_xticks([10,100])
ax[0][0].set_ylabel('b/a',fontsize=25)
ax[1][0].set_ylabel('c/a',fontsize=25)
ax[1][0].set_xlabel(r'Log[M$_{vir}$/M$_\odot$]',fontsize=25)
ax[1][1].set_xlabel(r'V$_{vir}$ [km/s]',fontsize=25)

Xsec,colors = ['CDM','SI3','SI10','vdXsec'],['k','r','b','g']
for i in [3]:
    Mvir = np.load(f'../DataFiles/Mvir.storm.{Xsec[i]}.z0.npy')
    Rvir = np.load(f'../DataFiles/Rvir.storm.{Xsec[i]}.z0.npy')
    Cont = np.load(f'../DataFiles/ContaminationFraction.storm.{Xsec[i]}.z0.npy')
    Shapes = pickle.load(open(f'../DataFiles/DarkMatterShapes.Top200.{Xsec[i]}.pickle','rb'))

    b,c,m,v,G = [],[],[],[],4.30091e-3 #pc Msol^-1 (km/s)^2

    for hnum in np.arange(200):
        if Cont[hnum]<.1:
            m.append(np.log10(Mvir[hnum]))
            v.append(np.sqrt((Mvir[hnum]*G)/(Rvir[hnum]*1000)))
            if not np.isnan(Shapes[str(hnum+1)]['rbins'][0]):
                i_inner = np.argmin(abs(Shapes[str(hnum+1)]['rbins'] - (Shapes[str(hnum+1)]['rbins'][-1]*.1 )))
                b.append(Shapes[str(hnum+1)]['b_pyn'][i_inner])
                c.append(Shapes[str(hnum+1)]['c_pyn'][i_inner])

    s = 1.5
    ax[0][0].scatter(m,b,s=s**2,c=colors[i],label=Xsec[i])
    ax[0][1].scatter(v,b,s=s**2,c=colors[i])
    ax[1][0].scatter(m,c,s=s**2,c=colors[i])
    ax[1][1].scatter(v,c,s=s**2,c=colors[i])

ax[0][0].legend(loc='lower right',prop={'size':20})
f.savefig('../Plots/ShapesVsEtc.png',bbox_inches='tight',pad_inches=.1)
meta = OSXMetaData('../Plots/ShapesVsEtc.png')
meta.creator='../Plots/ShapesVsEtc.Plots.py'