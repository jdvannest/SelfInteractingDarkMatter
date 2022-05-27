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
    M0 = np.load(f'../DataFiles/Mvir.N{args.npart}.CDM.{z}.npy')
    R0 = np.load(f'../DataFiles/Rvir.N{args.npart}.CDM.{z}.npy')
    M3 = np.load(f'../DataFiles/Mvir.N{args.npart}.SI3.{z}.npy')
    R3 = np.load(f'../DataFiles/Rvir.N{args.npart}.SI3.{z}.npy')
    M10 = np.load(f'../DataFiles/Mvir.N{args.npart}.SI10.{z}.npy')
    R10 = np.load(f'../DataFiles/Rvir.N{args.npart}.SI10.{z}.npy')
    if z in ['z0','z1']:
        MV = np.load(f'../DataFiles/Mvir.N{args.npart}.vdXsec.{z}.npy')
        RV = np.load(f'../DataFiles/Rvir.N{args.npart}.vdXsec.{z}.npy')
        MVc = np.load(f'../DataFiles/ContaminationFraction.N{args.npart}.vdXsec.{z}.npy')
    #Apply contamination fraction limit
    M0c = np.load(f'../DataFiles/ContaminationFraction.N{args.npart}.CDM.{z}.npy')
    M3c = np.load(f'../DataFiles/ContaminationFraction.N{args.npart}.SI3.{z}.npy')
    M10c = np.load(f'../DataFiles/ContaminationFraction.N{args.npart}.SI10.{z}.npy')
    contam_limit=.1
    M0 = M0[M0c<contam_limit]
    R0 = R0[M0c<contam_limit]
    M3 = M3[M3c<contam_limit]
    R3 = R3[M3c<contam_limit]
    M10 = M10[M10c<contam_limit]
    R10 = R10[M10c<contam_limit]
    if z in ['z0','z1']:
        MV = MV[MVc<contam_limit]
        RV = RV[MVc<contam_limit]
        pMV = np.zeros(len(mass_bins))
    pM0,pM3,pM10 = np.zeros(len(mass_bins)),np.zeros(len(mass_bins)),np.zeros(len(mass_bins))
    masses,profiles = [M0,M3,M10],[pM0,pM3,pM10]
    if z in ['z0','z1']:
        masses.append(MV)
        profiles.append(pMV)
    for i in np.arange(len(mass_bins)):
        for j in np.arange(len(masses)):
            for m in masses[j]:
                if np.log10(m)>mass_bins[i]: profiles[j][i]+=1
    
    ax[redshifts.index(z)][0].hist(np.log10(M0),mass_bins,histtype='step',color='k',label='CDM')
    ax[redshifts.index(z)][0].hist(np.log10(M3),mass_bins,histtype='step',color='r',label='SI3')
    ax[redshifts.index(z)][0].hist(np.log10(M10),mass_bins,histtype='step',color='b',label='SI10')

    ax[redshifts.index(z)][1].hist(np.sqrt((M0*G)/(R0*1000)),vel_bins,histtype='step',color='k',label='CDM')
    ax[redshifts.index(z)][1].hist(np.sqrt((M3*G)/(R3*1000)),vel_bins,histtype='step',color='r',label='SI3')
    ax[redshifts.index(z)][1].hist(np.sqrt((M10*G)/(R10*1000)),vel_bins,histtype='step',color='b',label='SI10')

    ax[redshifts.index(z)][2].plot(mass_bins,pM0/pM0[0],color='k',label='CDM')
    ax[redshifts.index(z)][2].plot(mass_bins,pM3/pM3[0],color='r',label='SI3')
    ax[redshifts.index(z)][2].plot(mass_bins,pM10/pM10[0],color='b',label='SI10')

    if z in ['z0','z1']:
        ax[redshifts.index(z)][0].hist(np.log10(MV),mass_bins,histtype='step',color='g',label='vdXsec')
        ax[redshifts.index(z)][1].hist(np.sqrt((MV*G)/(RV*1000)),vel_bins,histtype='step',color='g',label='vdXsec')
        ax[redshifts.index(z)][2].plot(mass_bins,pMV/pMV[0],color='g',label='vdXsec')
    
    ax[redshifts.index(z)][0].set_ylabel(f'{z}',fontsize=20)
    ax[redshifts.index(z)][2].set_ylabel(r'N(M$_{vir}>$M)',fontsize=20)

ax[4][0].set_xlabel('Log[M$_{vir}$/M$_\odot$]',fontsize=20)
ax[4][1].set_xlabel('V$_{vir}$ [km/s]',fontsize=20)
ax[4][2].set_xlabel('Log[M/M$_\odot$]',fontsize=20)
ax[0][2].legend(loc='lower left',prop={'size':15})

for r in np.arange(5):
    
    ax[r][0].set_xlim([5.5,11])
    ax[r][0].set_ylim([0.8,5000])
    ax[r][0].semilogy()
    ax[r][0].tick_params(axis='y',which='both',right=True,left=True)

    ax[r][1].set_xlim([.5,100])
    ax[r][1].semilogx()
    ax[r][1].set_ylim([0.8,5000])
    ax[r][1].semilogy()
    ax[r][1].tick_params(axis='y',which='both',labelleft=False,direction='in')

    ax[r][2].set_xlim([5.5,11])
    ax[r][2].set_ylim([2e-5,1.5])
    ax[r][2].semilogy()
    ax[r][2].yaxis.set_label_position("right")
    ax[r][2].yaxis.tick_right()

f.savefig(f'../Plots/RedshiftComparison.N{args.npart}.png',bbox_inches='tight',pad_inches=.1)