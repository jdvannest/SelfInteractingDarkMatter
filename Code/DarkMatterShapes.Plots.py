import argparse,pickle
import matplotlib.pylab as plt
from osxmetadata import OSXMetaData

parser = argparse.ArgumentParser()
parser.add_argument("-c","--cross_section",required=True,choices=['3','10','30','50'])
#parser.add_argument("-p","--pynbody",action="store_true")
args = parser.parse_args()

data = pickle.load(open(f'../DataFiles/DarkMatterShapes.SI{args.cross_section}.pickle','rb'))

b,c = [[],[]]
for halo in data:
    b.append(data[halo]['b'])
    c.append(data[halo]['c'])

f,ax = plt.subplots(1,1,figsize=(8,8))
ax.set_xlim([0,1])
ax.set_ylim([0,1])
ax.set_xlabel('b/a',fontsize=25)
ax.set_ylabel('c/a',fontsize=25)
ax.plot([0,1],[0,1],c='0.5',linestyle='--')
ax.scatter(b,c,c='k',marker='.',s=.75**2)
f.savefig(f'../Plots/DarkMatterShapes.SI{args.cross_section}.AHF.png',bbox_inches='tight',pad_inches=.1)
meta = OSXMetaData(f'../Plots/DarkMatterShapes.SI{args.cross_section}.AHF.png')
meta.creator='DarkMatterShapes.Plots.py'