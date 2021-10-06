import argparse,pickle
import numpy as np
import matplotlib.pylab as plt
import matplotlib.colors as colors
from osxmetadata import OSXMetaData
from scipy.stats import kde
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

b,c,bpyn,cpyn = [[],[],[],[]]
for cs in ['CDM','SI3','SI10']:
    b.append([])
    c.append([])
    bpyn.append([])
    cpyn.append([])
    filename = f'../DataFiles/DarkMatterShapes.{cs}.pickle'
    data = pickle.load(open(filename,'rb'))
    for halo in data:
        b[-1].append(data[halo]['b'])
        c[-1].append(data[halo]['c'])
        bpyn[-1].append(data[halo]['b_pyn'][-1])
        cpyn[-1].append(data[halo]['c_pyn'][-1])

bins = np.linspace(0,1,101)
f,ax = plt.subplots(2,3,figsize=(22,15))
plt.subplots_adjust(wspace=0,hspace=0)
for i in [0,1]:
    for j in [0,1,2]:
        ax[i][j].set_xticks([])
        ax[i][j].set_yticks([])
        ax[i][j].set_xlim([0,1])
        ax[i][j].set_ylim([0,1])
        ax[i][j].tick_params(length=5,labelsize=15)
        ax[i][j].plot([0,1],[0,1],c='0.5',linestyle='--')
        if i == 1: ax[i][j].set_xticks([0,.2,.4,.6,.8])
        if j == 0: ax[i][j].set_yticks([0,.2,.4,.6,.8])
ax[0][0].set_yticks([0,.2,.4,.6,.8,1])
ax[1][2].set_xticks([0,.2,.4,.6,.8,1])
ax[1][1].set_xlabel('b/a',fontsize=30)
ax[0][0].set_ylabel('c/a [AHF]',fontsize=30)
ax[1][0].set_ylabel('c/a [Pynbody]',fontsize=30)
ax[0][0].set_title('CDM',fontsize=25)
ax[0][1].set_title('SI3',fontsize=25)
ax[0][2].set_title('SI10',fontsize=25)

ax[0][0].hist2d(b[0],c[0],[bins,bins],cmap='Greys',density=True)
ax[0][1].hist2d(b[1],c[1],[bins,bins],cmap=truncate_colormap(plt.get_cmap('seismic'),0.5,1),density=True)
ax[0][2].hist2d(b[2],c[2],[bins,bins],cmap=truncate_colormap(plt.get_cmap('seismic_r'),0.5,1),density=True)
ax[1][0].hist2d(bpyn[0],cpyn[0],[bins,bins],cmap='Greys',density=True)
ax[1][1].hist2d(bpyn[1],cpyn[1],[bins,bins],cmap=truncate_colormap(plt.get_cmap('seismic'),0.5,1),density=True)
ax[1][2].hist2d(bpyn[2],cpyn[2],[bins,bins],cmap=truncate_colormap(plt.get_cmap('seismic_r'),0.5,1),density=True)

for i in [0,1,2]:
    ax2 = ax[0][i].twinx()
    ax2.set_ylim([0,50])
    ax2.set_yticks([])
    ax2.hist(b[0],bins,histtype='step',facecolor='None',edgecolor='k',density=True)
    ax3 = ax[0][i].twiny()
    ax3.set_xlim([0,50])
    ax3.set_xticks([])
    ax3.hist(c[0],bins,histtype='step',facecolor='None',edgecolor='k',density=True,orientation='horizontal')
    ax4 = ax[0][i].twinx()
    ax4.set_ylim([0,50])
    ax4.set_yticks([])
    ax4.hist(b[1],bins,histtype='step',facecolor='None',edgecolor='r',density=True)
    ax5 = ax[0][i].twiny()
    ax5.set_xlim([0,50])
    ax5.set_xticks([])
    ax5.hist(c[1],bins,histtype='step',facecolor='None',edgecolor='r',density=True,orientation='horizontal')
    ax6 = ax[0][i].twinx()
    ax6.set_ylim([0,50])
    ax6.set_yticks([])
    ax6.hist(b[2],bins,histtype='step',facecolor='None',edgecolor='b',density=True)
    ax7 = ax[0][i].twiny()
    ax7.set_xlim([0,50])
    ax7.set_xticks([])
    ax7.hist(c[2],bins,histtype='step',facecolor='None',edgecolor='b',density=True,orientation='horizontal')

for i in [0,1,2]:
    ax2 = ax[1][i].twinx()
    ax2.set_ylim([0,50])
    ax2.set_yticks([])
    ax2.hist(bpyn[0],bins,histtype='step',facecolor='None',edgecolor='k',density=True)
    ax3 = ax[1][i].twiny()
    ax3.set_xlim([0,50])
    ax3.set_xticks([])
    ax3.hist(cpyn[0],bins,histtype='step',facecolor='None',edgecolor='k',density=True,orientation='horizontal')
    ax4 = ax[1][i].twinx()
    ax4.set_ylim([0,50])
    ax4.set_yticks([])
    ax4.hist(bpyn[1],bins,histtype='step',facecolor='None',edgecolor='r',density=True)
    ax5 = ax[1][i].twiny()
    ax5.set_xlim([0,50])
    ax5.set_xticks([])
    ax5.hist(cpyn[1],bins,histtype='step',facecolor='None',edgecolor='r',density=True,orientation='horizontal')
    ax6 = ax[1][i].twinx()
    ax6.set_ylim([0,50])
    ax6.set_yticks([])
    ax6.hist(bpyn[2],bins,histtype='step',facecolor='None',edgecolor='b',density=True)
    ax7 = ax[1][i].twiny()
    ax7.set_xlim([0,50])
    ax7.set_xticks([])
    ax7.hist(cpyn[2],bins,histtype='step',facecolor='None',edgecolor='b',density=True,orientation='horizontal')

f.savefig(f'../Plots/DarkMatterShapes.MasterHistogram.png',bbox_inches='tight',pad_inches=.1)
meta = OSXMetaData(f'../Plots/DarkMatterShapes.MasterHistogram.png')
meta.creator='DarkMatterShapes.MasterHistogram.py'