import argparse
import numpy as np
import matplotlib.pylab as plt

parser = argparse.ArgumentParser()
parser.add_argument('-s','--simulation',choices=['storm','h148'],required=True)
parser.add_argument('-z','--redshift',choices=['z0','z1','z2','z3','z4','z10','z13.5'],required=True)
parser.add_argument('-n','--npart',default=64)
args = parser.parse_args()

#Load Data
N0 = np.load(f'../DataFiles/Npart.{args.simulation}.CDM.{args.redshift}.npy')
NV = np.load(f'../DataFiles/Npart.{args.simulation}.vdXsec.{args.redshift}.npy')
M0 = np.load(f'../DataFiles/Mvir.{args.simulation}.CDM.{args.redshift}.npy')
MV = np.load(f'../DataFiles/Mvir.{args.simulation}.vdXsec.{args.redshift}.npy')
C0 = np.load(f'../DataFiles/ContaminationFraction.{args.simulation}.CDM.{args.redshift}.npy')
CV = np.load(f'../DataFiles/ContaminationFraction.{args.simulation}.vdXsec.{args.redshift}.npy')
if args.simulation=='storm':
    N3 = np.load(f'../DataFiles/Npart.{args.simulation}.SI3.{args.redshift}.npy')
    N10 = np.load(f'../DataFiles/Npart.{args.simulation}.SI10.{args.redshift}.npy')
    M3 = np.load(f'../DataFiles/Mvir.{args.simulation}.SI3.{args.redshift}.npy')
    M10 = np.load(f'../DataFiles/Mvir.{args.simulation}.SI10.{args.redshift}.npy')
    C3 = np.load(f'../DataFiles/ContaminationFraction.{args.simulation}.SI3.{args.redshift}.npy')
    C10 = np.load(f'../DataFiles/ContaminationFraction.{args.simulation}.SI10.{args.redshift}.npy')
    if args.redshift in ['z2','z3','z4','z13.5']:
        N50 = np.load(f'../DataFiles/Npart.{args.simulation}.SI50.{args.redshift}.npy')
        M50 = np.load(f'../DataFiles/Mvir.{args.simulation}.SI50.{args.redshift}.npy')
        C50 = np.load(f'../DataFiles/ContaminationFraction.{args.simulation}.SI50.{args.redshift}.npy')
if args.redshift=='z13.5':
    NV = np.load(f'../DataFiles/Npart.{args.simulation}.VTS.{args.redshift}.npy')
    MV = np.load(f'../DataFiles/Mvir.{args.simulation}.VTS.{args.redshift}.npy')
    CV = np.load(f'../DataFiles/ContaminationFraction.{args.simulation}.VTS.{args.redshift}.npy')

#Apply particle count limit
M0 = M0[N0>(int(args.npart)-1)]
C0 = C0[N0>(int(args.npart)-1)]
MV = MV[NV>(int(args.npart)-1)]
CV = CV[NV>(int(args.npart)-1)]
if args.simulation=='storm':
    M3 = M3[N3>(int(args.npart)-1)]
    C3 = C3[N3>(int(args.npart)-1)]
    M10 = M10[N10>(int(args.npart)-1)]
    C10 = C10[N10>(int(args.npart)-1)]
    if args.redshift in ['z2','z3','z4','z13.5']:
        M50 = M50[N50>(int(args.npart)-1)]
        C50 = C50[N50>(int(args.npart)-1)]
#Apply contamination fraction limit
contam_limit=.1
M0 = M0[C0<contam_limit]
MV = MV[CV<contam_limit]
if args.simulation=='storm':
    M3 = M3[C3<contam_limit]
    M10 = M10[C10<contam_limit]
    if args.redshift in ['z2','z3','z4','z13.5']: M50 = M50[C50<contam_limit]

#mass_bins = np.linspace(6,13.5,100)
mass_bins = np.linspace(5.5,11,100)

if args.simulation=='storm':
    pM0,pM3,pM10,pMV = np.zeros(len(mass_bins)),np.zeros(len(mass_bins)),np.zeros(len(mass_bins)),np.zeros(len(mass_bins))
    iM0,iM3,iM10,iMV = np.zeros(len(mass_bins)),np.zeros(len(mass_bins)),np.zeros(len(mass_bins)),np.zeros(len(mass_bins))
    masses,profiles,inverted = [M0,M3,M10,MV],[pM0,pM3,pM10,pMV],[iM0,iM3,iM10,iMV]
    if args.redshift in ['z2','z3','z4','z13.5']:
        pM50,iM50 = np.zeros(len(mass_bins)),np.zeros(len(mass_bins))
        masses[-1] = M50
        profiles[-1] = pM50
        inverted[-1] = iM50
        masses.append(MV)
        profiles.append(pMV)
        inverted.append(iMV)
else:
    pM0,pMV = np.zeros(len(mass_bins)),np.zeros(len(mass_bins))
    iM0,iMV = np.zeros(len(mass_bins)),np.zeros(len(mass_bins))
    masses,profiles,inverted = [M0,MV],[pM0,pMV],[iM0,iMV]

for i in np.arange(len(mass_bins)):
    for j in np.arange(len(masses)):
        for m in masses[j]:
            if np.log10(m)>mass_bins[i]: profiles[j][i]+=1
            if np.log10(m)<mass_bins[i]: inverted[j][i]+=1

#Cumulative Mass Functions
for b in [False,True]:
    f,ax = plt.subplots(1,2,figsize=(12,5))
    ax[0].set_ylabel(r'N(M$_{vir}>$M)',fontsize=15)
    ax[1].set_ylabel(r'Normalized N(M$_{vir}>$M)',fontsize=15)
    if b:
        fname='.LogScale'
    else:
        fname=''
        ax[0].set_ylim([0,len(M0)])
        ax[1].set_ylim([0,1])
    for i in [0,1]:
        ax[i].set_xlabel(r'Log[M/M$_\odot$]',fontsize=15)
        ax[i].set_xlim([5.5,11])
        ax[i].tick_params(labelsize=10)
        if b: ax[i].semilogy()

    ax[0].plot(mass_bins,pM0,c='k',label='CDM')
    if args.simulation=='storm':
        ax[0].plot(mass_bins,pM3,c='r',label='SI3')
        ax[0].plot(mass_bins,pM10,c='b',label='SI10')
        if args.redshift in ['z2','z3','z4','z13.5']: ax[0].plot(mass_bins,pM50,c='orange',label='SI50')
    #ax[0].plot(mass_bins,pMV,c='g',label='vdXsec')

    ax[1].plot(mass_bins,pM0/pM0[0],c='k',label='CDM')
    if args.simulation=='storm':
        ax[1].plot(mass_bins,pM3/pM3[0],c='r',label='SI3')
        ax[1].plot(mass_bins,pM10/pM10[0],c='b',label='SI10')
        if args.redshift in ['z2','z3','z4','z13.5']: ax[1].plot(mass_bins,pM50/pM50[0],c='orange',label='SI50')
    #ax[1].plot(mass_bins,pMV/pMV[0],c='g',label='vdXsec')

    ax[0].legend(loc='upper right',prop={'size':12})

    f.savefig(f'../Plots/MassFunctionComparison.{args.simulation}.N{args.npart}.{args.redshift}{fname}.png',bbox_inches='tight',pad_inches=.1)

#Inverted Mass Functions
for b in [False,True]:
    f,ax = plt.subplots(1,2,figsize=(12,5))
    ax[0].set_ylabel(r'N(M$_{vir}<$M)',fontsize=15)
    ax[1].set_ylabel(r'Normalized N(M$_{vir}<$M)',fontsize=15)
    if b:
        fname='.LogScale'
    else:
        fname=''
        ax[0].set_ylim([0,len(M0)])
        ax[1].set_ylim([0,1])
    for i in [0,1]:
        ax[i].set_xlabel(r'Log[M/M$_\odot$]',fontsize=15)
        ax[i].set_xlim([5.5,11])
        ax[i].tick_params(labelsize=10)
        if b: ax[i].semilogy()

    ax[0].plot(mass_bins,iM0,c='k',label='CDM')
    if args.simulation=='storm':
        ax[0].plot(mass_bins,iM3,c='r',label='SI3')
        ax[0].plot(mass_bins,iM10,c='b',label='SI10')
        if args.redshift in ['z2','z3','z4','z13.5']: ax[0].plot(mass_bins,iM50,c='orange',label='SI50')
    #ax[0].plot(mass_bins,iMV,c='g',label='vdXsec')

    ax[1].plot(mass_bins,iM0/iM0[-1],c='k',label='CDM')
    if args.simulation=='storm':
        ax[1].plot(mass_bins,iM3/iM3[-1],c='r',label='SI3')
        ax[1].plot(mass_bins,iM10/iM10[-1],c='b',label='SI10')
        if args.redshift in ['z2','z3','z4','z13.5']: ax[1].plot(mass_bins,iM50/iM50[-1],c='orange',label='SI50')
    #ax[1].plot(mass_bins,iMV/iMV[-1],c='g',label='vdXsec')

    ax[0].legend(loc='upper left',prop={'size':12})

    f.savefig(f'../Plots/MassFunctionComparison.{args.simulation}.Inverted.N{args.npart}.{args.redshift}{fname}.png',bbox_inches='tight',pad_inches=.1)


#Mass Histograms
for b in [False,True]:
    f,ax = plt.subplots(1,1,figsize=(6,4))
    ax.set_xlabel(r'Log[M$_{vir}$/M$_\odot$]',fontsize=15)
    if b: 
        ax.set_ylabel('Normalized N',fontsize=15)
    else:
        ax.set_ylabel('N',fontsize=15)
    ax.tick_params(labelsize=10)
    ax.set_xlim([5.5,11])

    ax.hist(np.log10(M0),mass_bins,histtype='step',density=b,color='k',label='CDM')
    if args.simulation=='storm':
        ax.hist(np.log10(M3),mass_bins,histtype='step',density=b,color='r',label='SI3')
        ax.hist(np.log10(M10),mass_bins,histtype='step',density=b,color='b',label='SI10')
        if args.redshift in ['z2','z3','z4','z13.5']: ax.hist(np.log10(M50),mass_bins,histtype='step',density=b,color='orange',label='SI50')
    #ax.hist(np.log10(MV),mass_bins,histtype='step',density=b,color='g',label='vdXsec')

    if not b:
        ax.semilogy()
        fname = 'LogScale'
    else:
        fname = 'Normalized'
    ax.legend(loc='upper right',prop={'size':12})
    f.savefig(f'../Plots/MassHistogramComparison.{args.simulation}.N{args.npart}.{args.redshift}.{fname}.png',bbox_inches='tight',pad_inches=.1)