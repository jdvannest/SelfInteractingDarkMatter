import os,pickle
import numpy as np
import matplotlib.pylab as plt

os.system(f'python Config.py')
config = pickle.load(open('Config.pickle','rb'))

C,color,Zs,Hs = ['CDM','SI3','SI10','vdXsec'],['k','r','b','g'],[0,0,0,0],[0,0,0,0]

for c in C:
    path = config['storm']['z0']['simpaths'][c]
    os.chdir('/'.join(path.split('/')[:-1]))
    ld = os.listdir()
    tsteps = []
    pre = '.'.join(path.split('/')[-1].split('.')[:-1])+'.'
    for i in ld: 
        if i.startswith(pre) and len(i)==len(path.split('/')[-1]):
            tsteps.append(i.split('.')[-1])
    tsteps.sort()
    z,hnum,j = np.zeros(len(tsteps)),np.zeros(len(tsteps)),0
    for t in tsteps:
        for i in ld:
            if i.startswith(pre+t) and i.split('.')[-1]=='AHF_halos':
                z[j] = float(i.split('z')[1][:5])
                num = 0
                with open(i) as f:
                    L = f.readlines()
                    del L[0]
                for l in L:
                    if float(l.split('\t')[37])>.9: num+=1
                hnum[j] = num
        j+=1
    Zs[C.index(c)] = z
    Hs[C.index(c)] = hnum

f,ax=plt.subplots(1,1)
ax.set_xlabel('Log[Z]')
ax.set_ylabel(r'N$_{halo}$')
ax.semilogx()
ax.invert_xaxis()
for i in [0,1,2,3]:
    ax.plot(Zs[i],Hs[i],color=color[i],label=C[i])
ax.legend(loc='upper left')
f.savefig(f'/home1/08902/tg882017/work2/SelfInteractingDarkMatter/Plots/HaloCountHistory.png',bbox_inches='tight',pad_inches=.1)
ax.semilogy()
ax.set_ylabel(r'Log[N$_{halo}$]')
f.savefig(f'/home1/08902/tg882017/work2/SelfInteractingDarkMatter/Plots/HaloCountHistory.Log.png',bbox_inches='tight',pad_inches=.1)

zc = Zs[0]
Hc = Hs[0]
Hn = [np.zeros(len(zc)),np.zeros(len(zc)),np.zeros(len(zc))]
for j in [0,1,2]:
    for i in np.arange(len(zc)):
        if zc[i] in Zs[j+1]:
            Hn[j][i] = Hs[j+1][np.where(Zs[j+1]==zc[i])[0][0]]
        else:
            Hn[j][i] = (Hs[j+1][Zs[j+1]<zc[i]][-1] + Hs[j+1][Zs[j+1]>zc[i]][0])/2

f,ax=plt.subplots(1,1)
ax.set_xlabel('Log[Z]')
ax.set_ylabel(r'N$_{halo}$/N$_{halo,CDM}$')
ax.semilogx()
ax.invert_xaxis()
for i in [0,1,2]:
    ax.plot(zc,Hn[i]/Hc,color=color[i+1],label=C[i+1])
ax.legend(loc='upper left')
f.savefig(f'/home1/08902/tg882017/work2/SelfInteractingDarkMatter/Plots/HaloCountHistory.Normalized.png',bbox_inches='tight',pad_inches=.1)