import pickle
import numpy as np
import matplotlib.pylab as plt
from osxmetadata import OSXMetaData

f,ax = plt.subplots(2,2,figsize=(8,8))
plt.subplots_adjust(hspace=0,wspace=0)
for i in [0,1]:
    for j in [0,1]:
        ax[i][j].set_xlim([0,1])
        ax[i][j].set_ylim([0,1])
        ax[i][j].set_xticks([])
        ax[i][j].set_yticks([])
        ax[i][j].tick_params(length=5,labelsize=15)
ax[1][0].set_xticks([0,.2,.4,.6,.8])
ax[1][0].set_yticks([0,.2,.4,.6,.8])
ax[1][1].set_xticks([0,.2,.4,.6,.8,1])
ax[0][0].set_yticks([0,.2,.4,.6,.8,1])
ax[0][0].set_ylabel('b/a',fontsize=25)
ax[1][0].set_ylabel('c/a',fontsize=25)
ax[0][0].set_title(r'R$_{vir}$ [AHF]',fontsize=25)
ax[0][1].set_title(r'.1*R$_{vir}$ [Pyn]',fontsize=25)

Xsec,colors,bins = ['CDM','SI3','SI10','SI50','vdXsec'],['k','r','b','orange','g'],np.linspace(0,1,101)
for i in [0,1,2,3,4]:
    Pyn = pickle.load(open(f'../DataFiles/DarkMatterShapes.Top200.{Xsec[i]}.pickle','rb'))
    b,c,bi,ci = [],[],[],[]

    for halo in Pyn:
        if not np.isnan(Pyn[halo]['rbins'][0]) and int(halo)<201:
            b.append(Pyn[halo]['b_ahf'])
            c.append(Pyn[halo]['c_ahf'])
            i_inner = np.argmin(abs(Pyn[halo]['rbins'] - (Pyn[halo]['rbins'][-1]*.1 )))
            bi.append(Pyn[halo]['b_pyn'][i_inner])
            ci.append(Pyn[halo]['c_pyn'][i_inner])
    
    ax[0][0].hist(b,bins,histtype='step',density=True,cumulative=True,color=colors[i],label=Xsec[i])
    ax[0][1].hist(bi,bins,histtype='step',density=True,cumulative=True,color=colors[i])
    ax[1][0].hist(c,bins,histtype='step',density=True,cumulative=True,color=colors[i])
    ax[1][1].hist(ci,bins,histtype='step',density=True,cumulative=True,color=colors[i])

ax[0][0].legend(loc='upper left',prop={'size':20})
f.savefig('../Plots/DarkMatterShapes.Top.vdXsec.png',bbox_inches='tight',pad_inches=.1)
meta = OSXMetaData('../Plots/DarkMatterShapes.Top.vdXsec.png')
meta.creator='DarkMatterShapes.vdXsec.Plots.py'



#Pynbody for Rvir and .1*Rvir
f,ax = plt.subplots(2,2,figsize=(8,8))
plt.subplots_adjust(hspace=0,wspace=0)
for i in [0,1]:
    for j in [0,1]:
        ax[i][j].set_xlim([0,1])
        ax[i][j].set_ylim([0,1])
        ax[i][j].set_xticks([])
        ax[i][j].set_yticks([])
        ax[i][j].tick_params(length=5,labelsize=15)
ax[1][0].set_xticks([0,.2,.4,.6,.8])
ax[1][0].set_yticks([0,.2,.4,.6,.8])
ax[1][1].set_xticks([0,.2,.4,.6,.8,1])
ax[0][0].set_yticks([0,.2,.4,.6,.8,1])
ax[0][0].set_ylabel('b/a',fontsize=25)
ax[1][0].set_ylabel('c/a',fontsize=25)
ax[0][0].set_title(r'R$_{vir}$',fontsize=25)
ax[0][1].set_title(r'.1*R$_{vir}$',fontsize=25)

Xsec,colors,bins = ['CDM','SI3','SI10','SI50','vdXsec'],['k','r','b','orange','g'],np.linspace(0,1,101)
for i in [0,1,2,3,4]:
    Pyn = pickle.load(open(f'../DataFiles/DarkMatterShapes.Top200.{Xsec[i]}.pickle','rb'))
    b,c,bi,ci = [],[],[],[]

    for halo in Pyn:
        if not np.isnan(Pyn[halo]['rbins'][0]) and int(halo)<201:
            b.append(Pyn[halo]['b_pyn'][-1])
            c.append(Pyn[halo]['c_pyn'][-1])
            i_inner = np.argmin(abs(Pyn[halo]['rbins'] - (Pyn[halo]['rbins'][-1]*.1 )))
            bi.append(Pyn[halo]['b_pyn'][i_inner])
            ci.append(Pyn[halo]['c_pyn'][i_inner])
    
    ax[0][0].hist(b,bins,histtype='step',density=True,cumulative=True,color=colors[i],label=Xsec[i])
    ax[0][1].hist(bi,bins,histtype='step',density=True,cumulative=True,color=colors[i])
    ax[1][0].hist(c,bins,histtype='step',density=True,cumulative=True,color=colors[i])
    ax[1][1].hist(ci,bins,histtype='step',density=True,cumulative=True,color=colors[i])

ax[0][0].legend(loc='upper left',prop={'size':20})
f.savefig('../Plots/DarkMatterShapes.BothPynbody.Top.vdXsec.png',bbox_inches='tight',pad_inches=.1)
meta = OSXMetaData('../Plots/DarkMatterShapes.BothPynbody.Top.vdXsec.png')
meta.creator='DarkMatterShapes.vdXsec.Plots.py'