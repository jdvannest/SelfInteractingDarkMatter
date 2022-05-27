import argparse
import numpy as np
import matplotlib.pylab as plt

parser = argparse.ArgumentParser()
parser.add_argument('-z','--redshift',choices=['z0','z1','z2','z3','z4'],required=True)
parser.add_argument('-n','--npart',default=64)
args = parser.parse_args()

#Load Data
M0 = np.load(f'../DataFiles/Mvir.N{args.npart}.CDM.{args.redshift}.npy')
M3 = np.load(f'../DataFiles/Mvir.N{args.npart}.SI3.{args.redshift}.npy')
M10 = np.load(f'../DataFiles/Mvir.N{args.npart}.SI10.{args.redshift}.npy')
MV = np.load(f'../DataFiles/Mvir.N{args.npart}.vdXsec.{args.redshift}.npy')
#Apply contamination fraction limit
M0c = np.load(f'../DataFiles/ContaminationFraction.N{args.npart}.CDM.{args.redshift}.npy')
M3c = np.load(f'../DataFiles/ContaminationFraction.N{args.npart}.SI3.{args.redshift}.npy')
M10c = np.load(f'../DataFiles/ContaminationFraction.N{args.npart}.SI10.{args.redshift}.npy')
MVc = np.load(f'../DataFiles/ContaminationFraction.N{args.npart}.vdXsec.{args.redshift}.npy')
contam_limit=.1
M0 = M0[M0c<contam_limit]
M3 = M3[M3c<contam_limit]
M10 = M10[M10c<contam_limit]
MV = MV[MVc<contam_limit]

#mass_bins = np.linspace(6,13.5,100)
mass_bins = np.linspace(5.5,11,100)

pM0,pM3,pM10,pMV = np.zeros(len(mass_bins)),np.zeros(len(mass_bins)),np.zeros(len(mass_bins)),np.zeros(len(mass_bins))
iM0,iM3,iM10,iMV = np.zeros(len(mass_bins)),np.zeros(len(mass_bins)),np.zeros(len(mass_bins)),np.zeros(len(mass_bins))
masses,profiles,inverted = [M0,M3,M10,MV],[pM0,pM3,pM10,pMV],[iM0,iM3,iM10,iMV]

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
        ax[0].set_ylim([0,max([len(M0),len(M3),len(M10),len(MV)])])
        ax[1].set_ylim([0,1])
    for i in [0,1]:
        ax[i].set_xlabel(r'Log[M/M$_\odot$]',fontsize=15)
        ax[i].set_xlim([5.5,11])
        ax[i].tick_params(labelsize=10)
        if b: ax[i].semilogy()

    ax[0].plot(mass_bins,pM0,c='k',label='CDM')
    ax[0].plot(mass_bins,pM3,c='r',label='SI3')
    ax[0].plot(mass_bins,pM10,c='b',label='SI10')
    ax[0].plot(mass_bins,pMV,c='g',label='vdXsec')

    ax[1].plot(mass_bins,pM0/pM0[0],c='k',label='CDM')
    ax[1].plot(mass_bins,pM3/pM3[0],c='r',label='SI3')
    ax[1].plot(mass_bins,pM10/pM10[0],c='b',label='SI10')
    ax[1].plot(mass_bins,pMV/pMV[0],c='g',label='vdXsec')

    ax[0].legend(loc='upper right',prop={'size':12})

    f.savefig(f'../Plots/MassFunctionComparison.N{args.npart}.{args.redshift}{fname}.png',bbox_inches='tight',pad_inches=.1)

#Inverted Mass Functions
for b in [False,True]:
    f,ax = plt.subplots(1,2,figsize=(12,5))
    ax[0].set_ylabel(r'N(M$_{vir}<$M)',fontsize=15)
    ax[1].set_ylabel(r'Normalized N(M$_{vir}<$M)',fontsize=15)
    if b:
        fname='.LogScale'
    else:
        fname=''
        ax[0].set_ylim([0,max([len(M0),len(M3),len(M10),len(MV)])])
        ax[1].set_ylim([0,1])
    for i in [0,1]:
        ax[i].set_xlabel(r'Log[M/M$_\odot$]',fontsize=15)
        ax[i].set_xlim([5.5,11])
        ax[i].tick_params(labelsize=10)
        if b: ax[i].semilogy()

    ax[0].plot(mass_bins,iM0,c='k',label='CDM')
    ax[0].plot(mass_bins,iM3,c='r',label='SI3')
    ax[0].plot(mass_bins,iM10,c='b',label='SI10')
    ax[0].plot(mass_bins,iMV,c='g',label='vdXsec')

    ax[1].plot(mass_bins,iM0/iM0[-1],c='k',label='CDM')
    ax[1].plot(mass_bins,iM3/iM3[-1],c='r',label='SI3')
    ax[1].plot(mass_bins,iM10/iM10[-1],c='b',label='SI10')
    ax[1].plot(mass_bins,iMV/iMV[-1],c='g',label='vdXsec')

    ax[0].legend(loc='upper left',prop={'size':12})

    f.savefig(f'../Plots/MassFunctionComparison.Inverted.N{args.npart}.{args.redshift}{fname}.png',bbox_inches='tight',pad_inches=.1)


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
    ax.hist(np.log10(M3),mass_bins,histtype='step',density=b,color='r',label='SI3')
    ax.hist(np.log10(M10),mass_bins,histtype='step',density=b,color='b',label='SI10')
    ax.hist(np.log10(MV),mass_bins,histtype='step',density=b,color='g',label='vdXsec')

    if not b:
        ax.semilogy()
        fname = 'LogScale'
    else:
        fname = 'Normalized'
    ax.legend(loc='upper right',prop={'size':12})
    f.savefig(f'../Plots/MassHistogramComparison.N{args.npart}.{args.redshift}.{fname}.png',bbox_inches='tight',pad_inches=.1)