import argparse
import numpy as np
import matplotlib.pylab as plt

parser = argparse.ArgumentParser()
parser.add_argument('-n','--npart',default=64)
args = parser.parse_args()

redshifts = ['z0','z1','z2','z3','z4']
vel_bins = np.logspace(-.3,2.65,100)
mass_bins = np.linspace(5.5,11,100)
G = 4.30091e-3 #pc Msol^-1 (km/s)^2

f,ax = plt.subplots(5,3,figsize=(15,15))
plt.subplots_adjust(wspace=0,hspace=0)       
for z in redshifts:
    #Load Data
    N0 = np.load(f'../DataFiles/Npart.h148.CDM.{z}.npy')
    M0 = np.load(f'../DataFiles/Mvir.h148.CDM.{z}.npy')
    R0 = np.load(f'../DataFiles/Rvir.h148.CDM.{z}.npy')
    C0 = np.load(f'../DataFiles/ContaminationFraction.h148.CDM.{z}.npy')
    #Apply particle count limit
    M0 = M0[N0>(int(args.npart)-1)]
    C0 = C0[N0>(int(args.npart)-1)]
    R0 = R0[N0>(int(args.npart)-1)]
    #Apply contamination fraction limit
    contam_limit=.1
    M0 = M0[C0<contam_limit]
    R0 = R0[C0<contam_limit]
    if z in ['z1','z2','z3']:
        NV = np.load(f'../DataFiles/Npart.h148.vdXsec.{z}.npy')
        MV = np.load(f'../DataFiles/Mvir.h148.vdXsec.{z}.npy')
        RV = np.load(f'../DataFiles/Rvir.h148.vdXsec.{z}.npy')
        CV = np.load(f'../DataFiles/ContaminationFraction.h148.vdXsec.{z}.npy')
        MV = MV[NV>(int(args.npart)-1)]
        CV = CV[NV>(int(args.npart)-1)]
        RV = RV[NV>(int(args.npart)-1)]
        MV = MV[CV<contam_limit]
        RV = RV[CV<contam_limit]
    pM0,pMV = np.zeros(len(mass_bins)),np.zeros(len(mass_bins))
    masses,profiles = [M0],[pM0]
    if z in ['z1','z2','z3']:
        masses.append(MV)
        profiles.append(pMV)
    for i in np.arange(len(mass_bins)):
        for j in np.arange(len(masses)):
            for m in masses[j]:
                if np.log10(m)>mass_bins[i]: profiles[j][i]+=1
    
    ax[redshifts.index(z)][0].hist(np.log10(M0),mass_bins,histtype='step',color='k',label='CDM')
    ax[redshifts.index(z)][1].hist(np.sqrt((M0*G)/(R0*1000)),vel_bins,histtype='step',color='k',label='CDM')
    ax[redshifts.index(z)][2].plot(mass_bins,pM0,color='k',label='CDM')

    if z in ['z1','z2','z3']:
        ax[redshifts.index(z)][0].hist(np.log10(MV),mass_bins,histtype='step',color='g',label='vdXsec')
        ax[redshifts.index(z)][1].hist(np.sqrt((MV*G)/(RV*1000)),vel_bins,histtype='step',color='g',label='vdXsec')
        ax[redshifts.index(z)][2].plot(mass_bins,pMV,color='g',label='vdXsec')
    
    ax[redshifts.index(z)][0].set_ylabel(f'{z}',fontsize=20)
    ax[redshifts.index(z)][2].set_ylabel(r'N(M$_{vir}>$M)',fontsize=20)

ax[4][0].set_xlabel('Log[M$_{vir}$/M$_\odot$]',fontsize=20)
ax[4][1].set_xlabel('V$_{vir}$ [km/s]',fontsize=20)
ax[4][2].set_xlabel('Log[M/M$_\odot$]',fontsize=20)
ax[0][2].plot([],[],c='g',label='vdXsec')
ax[0][2].legend(loc='lower left',prop={'size':15})

for r in np.arange(5):
    
    ax[r][0].set_xlim([5.5,11])
    ax[r][0].set_ylim([0.8,500])
    ax[r][0].semilogy()
    ax[r][0].tick_params(axis='y',which='both',right=True,left=True)

    ax[r][1].set_xlim([.5,100])
    ax[r][1].semilogx()
    ax[r][1].set_ylim([0.8,500])
    ax[r][1].semilogy()
    ax[r][1].tick_params(axis='y',which='both',labelleft=False,direction='in')

    ax[r][2].set_xlim([5.5,11])
    ax[r][2].set_ylim([.8,1e4])
    ax[r][2].semilogy()
    ax[r][2].yaxis.set_label_position("right")
    ax[r][2].yaxis.tick_right()

f.savefig(f'../Plots/RedshiftComparison.h148.N{args.npart}.png',bbox_inches='tight',pad_inches=.1)