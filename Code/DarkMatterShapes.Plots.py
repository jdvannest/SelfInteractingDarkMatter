import argparse,pickle
import numpy as np
import matplotlib.pylab as plt
from osxmetadata import OSXMetaData
from scipy.stats import kde

parser = argparse.ArgumentParser()
parser.add_argument("-c","--cross_section",required=True,choices=['cdm','3','10','30','50'])
#parser.add_argument("-p","--pynbody",action="store_true")
args = parser.parse_args()

if args.cross_section == 'cdm':
    filename = f'../DataFiles/DarkMatterShapes.CDM.pickle'
    im = 'CDM'
else:
    filename = f'../DataFiles/DarkMatterShapes.SI{args.cross_section}.pickle'
    im = f'SI{args.cross_section}'

data = pickle.load(open(filename,'rb'))

b,c = [[],[]]
for halo in data:
    b.append(data[halo]['b'])
    c.append(data[halo]['c'])

f,ax = plt.subplots(1,1,figsize=(8,8))
ax.set_xlim([0,1])
ax.set_ylim([0,1])
ax.tick_params(length=5,labelsize=15)
ax.set_xlabel('b/a',fontsize=25)
ax.set_ylabel('c/a',fontsize=25)
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
ax.set_xlabel('b/a',fontsize=25)
ax.set_ylabel('c/a',fontsize=25)
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
ax.set_xlabel('b/a',fontsize=25)
ax.set_ylabel('c/a',fontsize=25)
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