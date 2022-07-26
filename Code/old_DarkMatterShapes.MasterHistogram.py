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
parser = argparse.ArgumentParser()
parser.add_argument("-t","--top",action="store_true")
parser.add_argument("-i","--inner",action="store_true")
args = parser.parse_args()

topfile = '.Top200' if args.top else ''
topapp = '.top' if args.top else ''
app = '.Inner' if args.inner else ''

b,c,bpyn,cpyn = [[],[],[],[]]
for cs in ['CDM','SI3','SI10']:
    b.append([])
    c.append([])
    bpyn.append([])
    cpyn.append([])
    filename = f'../DataFiles/DarkMatterShapes{topfile}.{cs}.pickle'
    data = pickle.load(open(filename,'rb'))
    halo_list = [h for h in data]
    for halo in halo_list:
        b[-1].append(data[halo]['b'])
        c[-1].append(data[halo]['c'])
        if not np.isnan(data[halo]['rbins'][0]):
            if args.inner:
                i_inner = np.argmin(abs(data[halo]['rbins'] - (data[halo]['rbins'][-1]*.1 )))
                bpyn[-1].append(data[halo]['b_pyn'][i_inner])
                cpyn[-1].append(data[halo]['c_pyn'][i_inner])
            else:
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

f.savefig(f'../Plots/DarkMatterShapes{topapp}{app}.MasterHistogram.png',bbox_inches='tight',pad_inches=.1)
meta = OSXMetaData(f'../Plots/DarkMatterShapes{topapp}{app}.MasterHistogram.png')
meta.creator='DarkMatterShapes.MasterHistogram.py'


f,ax = plt.subplots(2,2,figsize=(10,10))
plt.subplots_adjust(hspace=0,wspace=0)
for i in [0,1]:
    for j in [0,1]:
        ax[i][j].set_xlim([0,1])
        ax[i][j].set_ylim([0,1])
        ax[i][j].tick_params(length=5,labelsize=15)
        if j==1: ax[i][j].set_yticks([])
        if i==0: ax[i][j].set_xticks([])
ax[1][0].set_xticks([0,.2,.4,.6,.8])
ax[1][0].set_yticks([0,.2,.4,.6,.8])
ax[1][1].set_xticks([0,.2,.4,.6,.8,1])
ax[0][0].set_yticks([0,.2,.4,.6,.8,1])
ax[0][0].set_ylabel('b/a',fontsize=25)
ax[1][0].set_ylabel('c/a',fontsize=25)
ax[0][0].set_title('AHF',fontsize=25)
ax[0][1].set_title('Pynbody',fontsize=25)

ax[0][0].hist(b[0],bins,facecolor='None',edgecolor='k',histtype='step',density=True,cumulative=True,label='CDM')
ax[0][0].hist(b[1],bins,facecolor='None',edgecolor='r',histtype='step',density=True,cumulative=True,label='SI3')
ax[0][0].hist(b[2],bins,facecolor='None',edgecolor='b',histtype='step',density=True,cumulative=True,label='SI10')

ax[1][0].hist(c[0],bins,facecolor='None',edgecolor='k',histtype='step',density=True,cumulative=True)
ax[1][0].hist(c[1],bins,facecolor='None',edgecolor='r',histtype='step',density=True,cumulative=True)
ax[1][0].hist(c[2],bins,facecolor='None',edgecolor='b',histtype='step',density=True,cumulative=True)

ax[0][1].hist(bpyn[0],bins,facecolor='None',edgecolor='k',histtype='step',density=True,cumulative=True)
ax[0][1].hist(bpyn[1],bins,facecolor='None',edgecolor='r',histtype='step',density=True,cumulative=True)
ax[0][1].hist(bpyn[2],bins,facecolor='None',edgecolor='b',histtype='step',density=True,cumulative=True)

ax[1][1].hist(cpyn[0],bins,facecolor='None',edgecolor='k',histtype='step',density=True,cumulative=True)
ax[1][1].hist(cpyn[1],bins,facecolor='None',edgecolor='r',histtype='step',density=True,cumulative=True)
ax[1][1].hist(cpyn[2],bins,facecolor='None',edgecolor='b',histtype='step',density=True,cumulative=True)

ax[0][0].legend(loc='upper left',prop={'size':15})
f.savefig(f'../Plots/DarkMatterShapes{topapp}{app}.CDF.png',bbox_inches='tight',pad_inches=.1)
meta = OSXMetaData(f'../Plots/DarkMatterShapes{topapp}{app}.CDF.png')
meta.creator='DarkMatterShapes.MasterHistogram.py'



f,ax = plt.subplots(1,3,figsize=(22,7))
plt.subplots_adjust(wspace=0)
pre='Inner ' if args.inner else ''
for j in [0,1,2]:
    ax[j].set_xticks([0,.2,.4,.6,.8])
    ax[j].set_yticks([])
    ax[j].set_xlim([0,1])
    ax[j].set_ylim([0,1])
    ax[j].tick_params(length=5,labelsize=20)
    ax[j].plot([0,1],[0,1],c='0.5',linestyle='--')
ax[0].set_yticks([0,.2,.4,.6,.8,1])
ax[2].set_xticks([0,.2,.4,.6,.8,1])
ax[0].set_ylabel(f'{pre}c/a',fontsize=30)
ax[1].set_xlabel(f'{pre}b/a',fontsize=30)
ax[0].set_title('CDM',fontsize=25)
ax[1].set_title('SI3',fontsize=25)
ax[2].set_title('SI10',fontsize=25)

bins = np.linspace(0,1,101)
ax[0].hist2d(bpyn[0],cpyn[0],[bins,bins],cmap='Greys',density=True)
ax[1].hist2d(bpyn[1],cpyn[1],[bins,bins],cmap=truncate_colormap(plt.get_cmap('seismic'),0.5,1),density=True)
ax[2].hist2d(bpyn[2],cpyn[2],[bins,bins],cmap=truncate_colormap(plt.get_cmap('seismic_r'),0.5,1),density=True)

f.savefig(f'/Users/jdvannest/Desktop/PynbodyOnly.DarkMatterShapes{topapp}{app}.MasterHistogram.png',bbox_inches='tight',pad_inches=.1)
meta = OSXMetaData(f'/Users/jdvannest/Desktop/PynbodyOnly.DarkMatterShapes{topapp}{app}.MasterHistogram.png')
meta.creator='DarkMatterShapes.MasterHistogram.py'


f,ax = plt.subplots(1,2,figsize=(10,5))
plt.subplots_adjust(wspace=0)
for i in [0,1]:
    ax[i].set_xlim([0,1])
    ax[i].set_ylim([0,1])
    ax[i].set_yticks([])
    ax[i].set_xticks([0,.2,.4,.6,.8])
    ax[i].tick_params(length=5,labelsize=15)
ax[0].set_yticks([0,.2,.4,.6,.8,1])
ax[1].set_xticks([0,.2,.4,.6,.8,1])
ax[0].set_title(f'{pre}b/a',fontsize=25)
ax[1].set_title(f'{pre}c/a',fontsize=25)

ax[0].hist(bpyn[0],bins,facecolor='None',edgecolor='k',histtype='step',density=True,cumulative=True,label='CDM')
ax[0].hist(bpyn[1],bins,facecolor='None',edgecolor='r',histtype='step',density=True,cumulative=True,label='SI3')
ax[0].hist(bpyn[2],bins,facecolor='None',edgecolor='b',histtype='step',density=True,cumulative=True,label='SI10')

ax[1].hist(cpyn[0],bins,facecolor='None',edgecolor='k',histtype='step',density=True,cumulative=True)
ax[1].hist(cpyn[1],bins,facecolor='None',edgecolor='r',histtype='step',density=True,cumulative=True)
ax[1].hist(cpyn[2],bins,facecolor='None',edgecolor='b',histtype='step',density=True,cumulative=True)

ax[0].legend(loc='upper left',prop={'size':15})

f.savefig(f'/Users/jdvannest/Desktop/PynbodyOnly.DarkMatterShapes{topapp}{app}.CDF.png',bbox_inches='tight',pad_inches=.1)
meta = OSXMetaData(f'/Users/jdvannest/Desktop/PynbodyOnly.DarkMatterShapes{topapp}{app}.CDF.png')
meta.creator='DarkMatterShapes.MasterHistogram.py'