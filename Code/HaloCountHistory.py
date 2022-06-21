import os,pickle
import matplotlib.pylab as plt

os.system(f'python Config.py')
config = pickle.load(open('Config.pickle','rb'))

f,ax=plt.subplots(1,1)
ax.set_xlabel('Log[Z]')
ax.set_ylabel(r'N$_{halo}$')
ax.semilogx()
ax.invert_xaxis()

C,color = ['CDM','SI3','SI10','vdXsec'],['k','r','b','g']

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
    z,hnum = [],[]
    for t in tsteps:
        for i in ld:
            if i.startswith(pre+t) and i.split('.')[-1]=='AHF_halos':
                z.append(float(i.split('z')[1][:5]))
                hnum.append(int(os.popen(f'wc -l {i}').read().split(' ')[0]))
    ax.plot(z,hnum,color=color[C.index(c)],label=c)

ax.legend(loc='upper left')
f.savefig(f'/home1/08902/tg882017/work2/SelfInteractingDarkMatter/Plots/HaloCountHistory.png',bbox_inches='tight',pad_inches=.1)
ax.semilogy()
ax.set_ylabel(r'Log[N$_{halo}$]')
f.savefig(f'/home1/08902/tg882017/work2/SelfInteractingDarkMatter/Plots/HaloCountHistory.Log.png',bbox_inches='tight',pad_inches=.1)