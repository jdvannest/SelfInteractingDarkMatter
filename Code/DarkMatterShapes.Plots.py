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
parser.add_argument("-c","--cross_section",required=True,choices=['cdm','3','10','30','50'])
parser.add_argument("-p","--pynbody",action="store_true")
args = parser.parse_args()

if args.cross_section == 'cdm':
    filename = f'../DataFiles/DarkMatterShapes.CDM.pickle'
    im = 'CDM'
else:
    filename = f'../DataFiles/DarkMatterShapes.SI{args.cross_section}.pickle'
    im = f'SI{args.cross_section}'

data = pickle.load(open(filename,'rb'))

b,c,bpyn,cpyn = [[],[],[],[]]
for halo in data:
    b.append(data[halo]['b'])
    c.append(data[halo]['c'])
    bpyn.append(data[halo]['b_pyn'][-1])
    cpyn.append(data[halo]['c_pyn'][-1])

f,ax = plt.subplots(1,1,figsize=(8,8))
ax.set_xlim([0,1])
ax.set_ylim([0,1])
ax.tick_params(length=5,labelsize=15)
ax.set_xlabel('b/a [AHF]',fontsize=25)
ax.set_ylabel('c/a [AHF]',fontsize=25)
ax.set_title(im,fontsize=25)
ax.plot([0,1],[0,1],c='0.5',linestyle='--')
ax.scatter(b,c,c='k',marker='.',s=.75**2)
f.savefig(f'../Plots/DarkMatterShapes.{im}.AHF.png',bbox_inches='tight',pad_inches=.1)
meta = OSXMetaData(f'../Plots/DarkMatterShapes.{im}.AHF.png')
meta.creator='DarkMatterShapes.Plots.py'

nbins=100
bins = np.linspace(0,1,101)
k = kde.gaussian_kde([b,c])
xk, yk = np.mgrid[min(b):max(b):nbins*1j, min(c):max(c):nbins*1j]
zk = k(np.vstack([xk.flatten(), yk.flatten()]))

f,ax = plt.subplots(1,1,figsize=(8,8))
ax.set_xlim([0,1])
ax.set_ylim([0,1])
ax.tick_params(length=5,labelsize=15)
ax.set_xlabel('b/a [AHF]',fontsize=25)
ax.set_ylabel('c/a [AHF]',fontsize=25)
ax.set_title(im,fontsize=25)
ax.plot([0,1],[0,1],c='0.5',linestyle='--')
ax.scatter(b,c,c='k',marker='.',s=.75**2)
ax.contour(xk,yk,zk.reshape(xk.shape),levels=3,cmap='spring')
f.savefig(f'../Plots/DarkMatterShapes.{im}.AHF.Contour.png',bbox_inches='tight',pad_inches=.1)
meta = OSXMetaData(f'../Plots/DarkMatterShapes.{im}.AHF.Contour.png')
meta.creator='DarkMatterShapes.Plots.py'


f,ax = plt.subplots(1,1,figsize=(8,8))
ax.set_xlim([0,1])
ax.set_ylim([0,1])
ax.tick_params(length=5,labelsize=15)
ax.set_xlabel('b/a [AHF]',fontsize=25)
ax.set_ylabel('c/a [AHF]',fontsize=25)
ax.set_title(im,fontsize=25)
ax.plot([0,1],[0,1],c='0.5',linestyle='--')
ax.hist2d(b,c,[bins,bins],cmap='Greys',density=True)
ax2 = ax.twinx()
ax2.set_ylim([0,50])
ax2.set_yticks([])
ax2.hist(b,bins,histtype='step',facecolor='None',edgecolor='k',density=True)
ax3 = ax.twiny()
ax3.set_xlim([0,50])
ax3.set_xticks([])
ax3.hist(c,bins,histtype='step',facecolor='None',edgecolor='k',density=True,orientation='horizontal')
f.savefig(f'../Plots/DarkMatterShapes.{im}.AHF.Histogram.png',bbox_inches='tight',pad_inches=.1)
meta = OSXMetaData(f'../Plots/DarkMatterShapes.{im}.AHF.Histogram.png')
meta.creator='DarkMatterShapes.Plots.py'

if args.pynbody:
    f,ax = plt.subplots(1,1,figsize=(8,8))
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.tick_params(length=5,labelsize=15)
    ax.set_xlabel('b/a [Pynbody]',fontsize=25)
    ax.set_ylabel('c/a [Pynbody]',fontsize=25)
    ax.set_title(im,fontsize=25)
    ax.plot([0,1],[0,1],c='0.5',linestyle='--')
    ax.scatter(bpyn,cpyn,c='k',marker='.',s=.75**2)
    f.savefig(f'../Plots/DarkMatterShapes.{im}.Pynbody.png',bbox_inches='tight',pad_inches=.1)
    meta = OSXMetaData(f'../Plots/DarkMatterShapes.{im}.Pynbody.png')
    meta.creator='DarkMatterShapes.Plots.py'

    nbins=100
    bins = np.linspace(0,1,101)
    k = kde.gaussian_kde([bpyn,cpyn])
    xk, yk = np.mgrid[min(bpyn):max(bpyn):nbins*1j, min(cpyn):max(cpyn):nbins*1j]
    zk = k(np.vstack([xk.flatten(), yk.flatten()]))

    f,ax = plt.subplots(1,1,figsize=(8,8))
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.tick_params(length=5,labelsize=15)
    ax.set_xlabel('b/a [Pynbody]',fontsize=25)
    ax.set_ylabel('c/a [Pynbody]',fontsize=25)
    ax.set_title(im,fontsize=25)
    ax.plot([0,1],[0,1],c='0.5',linestyle='--')
    ax.scatter(bpyn,cpyn,c='k',marker='.',s=.75**2)
    ax.contour(xk,yk,zk.reshape(xk.shape),levels=3,cmap='spring')
    f.savefig(f'../Plots/DarkMatterShapes.{im}.Pynbody.Contour.png',bbox_inches='tight',pad_inches=.1)
    meta = OSXMetaData(f'../Plots/DarkMatterShapes.{im}.Pynbody.Contour.png')
    meta.creator='DarkMatterShapes.Plots.py'


    f,ax = plt.subplots(1,1,figsize=(8,8))
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.tick_params(length=5,labelsize=15)
    ax.set_xlabel('b/a [Pynbody]',fontsize=25)
    ax.set_ylabel('c/a [Pynbody]',fontsize=25)
    ax.set_title(im,fontsize=25)
    ax.plot([0,1],[0,1],c='0.5',linestyle='--')
    ax.hist2d(bpyn,cpyn,[bins,bins],cmap='Greys',density=True)
    ax2 = ax.twinx()
    ax2.set_ylim([0,50])
    ax2.set_yticks([])
    ax2.hist(bpyn,bins,histtype='step',facecolor='None',edgecolor='k',density=True)
    ax3 = ax.twiny()
    ax3.set_xlim([0,50])
    ax3.set_xticks([])
    ax3.hist(cpyn,bins,histtype='step',facecolor='None',edgecolor='k',density=True,orientation='horizontal')
    f.savefig(f'../Plots/DarkMatterShapes.{im}.Pynbody.Histogram.png',bbox_inches='tight',pad_inches=.1)
    meta = OSXMetaData(f'../Plots/DarkMatterShapes.{im}.Pynbody.Histogram.png')
    meta.creator='DarkMatterShapes.Plots.py'

    f,ax = plt.subplots(1,2,figsize=(15,8))
    plt.subplots_adjust(wspace=0)
    for i in [0,1]:
        ax[i].set_xlim([0,1])
        ax[i].set_ylim([0,1])
        ax[i].tick_params(length=5,labelsize=15)
        ax[i].plot([0,1],[0,1],c='0.5',linestyle='--')
    ax[0].set_xlabel('b/a [AHF]',fontsize=25)
    ax[0].set_ylabel('c/a [AHF]',fontsize=25)
    ax[1].set_xlabel('b/a [Pynbody]',fontsize=25)
    ax[1].set_ylabel('c/a [Pynbody]',fontsize=25)
    ax[0].set_xticks([0,.2,.4,.6,.8])
    ax[1].set_yticks([])
    ax[1].yaxis.set_label_position("right")
    f.suptitle(im,fontsize=25)
    ax[0].hist2d(b,c,[bins,bins],cmap='Greys',density=True)
    ax[1].hist2d(bpyn,cpyn,[bins,bins],cmap=truncate_colormap(plt.get_cmap('seismic'), 0.5, 1),density=True)
    ax02 = ax[0].twinx()
    ax02.set_ylim([0,50])
    ax02.set_yticks([])
    ax02.hist(b,bins,histtype='step',facecolor='None',edgecolor='k',density=True)
    ax03 = ax[0].twiny()
    ax03.set_xlim([0,50])
    ax03.set_xticks([])
    ax03.hist(c,bins,histtype='step',facecolor='None',edgecolor='k',density=True,orientation='horizontal')
    ax04 = ax[0].twinx()
    ax04.set_ylim([0,50])
    ax04.set_yticks([])
    ax04.hist(bpyn,bins,histtype='step',facecolor='None',edgecolor='r',density=True)
    ax05 = ax[0].twiny()
    ax05.set_xlim([0,50])
    ax05.set_xticks([])
    ax05.hist(cpyn,bins,histtype='step',facecolor='None',edgecolor='r',density=True,orientation='horizontal')
    ax12 = ax[1].twinx()
    ax12.set_ylim([0,50])
    ax12.set_yticks([])
    ax12.hist(b,bins,histtype='step',facecolor='None',edgecolor='k',density=True)
    ax13 = ax[1].twiny()
    ax13.set_xlim([0,50])
    ax13.set_xticks([])
    ax13.hist(c,bins,histtype='step',facecolor='None',edgecolor='k',density=True,orientation='horizontal')
    ax14 = ax[1].twinx()
    ax14.set_ylim([0,50])
    ax14.set_yticks([])
    ax14.hist(bpyn,bins,histtype='step',facecolor='None',edgecolor='r',density=True)
    ax15 = ax[1].twiny()
    ax15.set_xlim([0,50])
    ax15.set_xticks([])
    ax15.hist(cpyn,bins,histtype='step',facecolor='None',edgecolor='r',density=True,orientation='horizontal')
    f.savefig(f'../Plots/DarkMatterShapes.{im}.Comparison.Histogram.png',bbox_inches='tight',pad_inches=.1)
    meta = OSXMetaData(f'../Plots/DarkMatterShapes.{im}.Comparison.Histogram.png')
    meta.creator='DarkMatterShapes.Plots.py'