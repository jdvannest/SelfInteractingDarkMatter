import numpy as np
import matplotlib.pylab as plt

M0 = np.load('../DataFiles/masses.CDM.Mine.npy')
M3 = np.load('../DataFiles/masses.SI3.Mine.npy')
M10 = np.load('../DataFiles/masses.SI10.Mine.npy')
MV = np.load('../DataFiles/masses.vdXsec.Mine.npy')
O0 = np.load('../DataFiles/masses.CDM.Old.npy')
O3 = np.load('../DataFiles/masses.SI3.Old.npy')
O10 = np.load('../DataFiles/masses.SI10.Old.npy')

mass_bins = np.linspace(6,13.5,100)

pM0,pM3,pM10,pMV,pO0,pO3,pO10 = [np.zeros(len(mass_bins)),np.zeros(len(mass_bins)),np.zeros(len(mass_bins)),
                                np.zeros(len(mass_bins)),np.zeros(len(mass_bins)),np.zeros(len(mass_bins)),np.zeros(len(mass_bins))]
iM0,iM3,iM10,iMV,iO0,iO3,iO10 = [np.zeros(len(mass_bins)),np.zeros(len(mass_bins)),np.zeros(len(mass_bins)),
                                np.zeros(len(mass_bins)),np.zeros(len(mass_bins)),np.zeros(len(mass_bins)),np.zeros(len(mass_bins))]
masses,profiles,inverted = [M0,M3,M10,MV,O0,O3,O10],[pM0,pM3,pM10,pMV,pO0,pO3,pO10],[iM0,iM3,iM10,iMV,iO0,iO3,iO10]

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
        ax[0].set_ylim([0,10000])
        ax[1].set_ylim([0,1])
    for i in [0,1]:
        ax[i].set_xlabel(r'Log[M/M$_\odot$]',fontsize=15)
        ax[i].set_xlim([6,13.5])
        ax[i].tick_params(labelsize=10)
        if b: ax[i].semilogy()
    ax[0].plot(0,0,c='0.5',label='My AHF')
    ax[0].plot(0,0,c='0.5',linestyle='--',label='Old AHF')

    ax[0].plot(mass_bins,pM0,c='k',label='CDM')
    ax[0].plot(mass_bins,pO0,c='k',linestyle='--')
    ax[0].plot(mass_bins,pM3,c='r',label='SI3')
    ax[0].plot(mass_bins,pO3,c='r',linestyle='--')
    ax[0].plot(mass_bins,pM10,c='b',label='SI10')
    ax[0].plot(mass_bins,pO10,c='b',linestyle='--')
    ax[0].plot(mass_bins,pMV,c='g',label='vdXsec')

    ax[1].plot(mass_bins,pM0/pM0[0],c='k',label='CDM')
    ax[1].plot(mass_bins,pO0/pO0[0],c='k',linestyle='--')
    ax[1].plot(mass_bins,pM3/pM3[0],c='r',label='SI3')
    ax[1].plot(mass_bins,pO3/pO3[0],c='r',linestyle='--')
    ax[1].plot(mass_bins,pM10/pM10[0],c='b',label='SI10')
    ax[1].plot(mass_bins,pO10/pO10[0],c='b',linestyle='--')
    ax[1].plot(mass_bins,pMV/pMV[0],c='g',label='vdXsec')

    ax[0].legend(loc='upper right',prop={'size':12})

    f.savefig(f'../Plots/MassFunctionComparison{fname}.png',bbox_inches='tight',pad_inches=.1)

#Inverted Mass Functions
for b in [False,True]:
    f,ax = plt.subplots(1,2,figsize=(12,5))
    ax[0].set_ylabel(r'N(M$_{vir}<$M)',fontsize=15)
    ax[1].set_ylabel(r'Normalized N(M$_{vir}<$M)',fontsize=15)
    if b:
        fname='.LogScale'
    else:
        fname=''
        ax[0].set_ylim([0,10000])
        ax[1].set_ylim([0,1])
    for i in [0,1]:
        ax[i].set_xlabel(r'Log[M/M$_\odot$]',fontsize=15)
        ax[i].set_xlim([6,13.5])
        ax[i].tick_params(labelsize=10)
        if b: ax[i].semilogy()
    ax[0].plot(0,0,c='0.5',label='My AHF')
    ax[0].plot(0,0,c='0.5',linestyle='--',label='Old AHF')

    ax[0].plot(mass_bins,iM0,c='k',label='CDM')
    ax[0].plot(mass_bins,iO0,c='k',linestyle='--')
    ax[0].plot(mass_bins,iM3,c='r',label='SI3')
    ax[0].plot(mass_bins,iO3,c='r',linestyle='--')
    ax[0].plot(mass_bins,iM10,c='b',label='SI10')
    ax[0].plot(mass_bins,iO10,c='b',linestyle='--')
    ax[0].plot(mass_bins,iMV,c='g',label='vdXsec')

    ax[1].plot(mass_bins,iM0/iM0[-1],c='k',label='CDM')
    ax[1].plot(mass_bins,iO0/iO0[-1],c='k',linestyle='--')
    ax[1].plot(mass_bins,iM3/iM3[-1],c='r',label='SI3')
    ax[1].plot(mass_bins,iO3/iO3[-1],c='r',linestyle='--')
    ax[1].plot(mass_bins,iM10/iM10[-1],c='b',label='SI10')
    ax[1].plot(mass_bins,iO10/iO10[-1],c='b',linestyle='--')
    ax[1].plot(mass_bins,iMV/iMV[-1],c='g',label='vdXsec')

    ax[0].legend(loc='upper left',prop={'size':12})

    f.savefig(f'../Plots/MassFunctionComparison.Inverted{fname}.png',bbox_inches='tight',pad_inches=.1)


#Mass Histograms
for b in [False,True]:
    f,ax = plt.subplots(1,1,figsize=(6,4))
    ax.set_xlabel(r'Log[M$_{vir}$/M$_\odot$]',fontsize=15)
    if b: 
        ax.set_ylabel('Normalized N',fontsize=15)
    else:
        ax.set_ylabel('N',fontsize=15)
    ax.tick_params(labelsize=10)
    ax.set_xlim([6,13.5])
    ax.plot(0,0,c='0.5',label='My AHF')
    ax.plot(0,0,c='0.5',linestyle='--',label='Old AHF')

    ax.hist(np.log10(M0),mass_bins,histtype='step',density=b,color='k',label='CDM')
    ax.hist(np.log10(O0),mass_bins,histtype='step',density=b,color='k',linestyle='--')
    ax.hist(np.log10(M3),mass_bins,histtype='step',density=b,color='r',label='SI3')
    ax.hist(np.log10(O3),mass_bins,histtype='step',density=b,color='r',linestyle='--')
    ax.hist(np.log10(M10),mass_bins,histtype='step',density=b,color='b',label='SI10')
    ax.hist(np.log10(O10),mass_bins,histtype='step',density=b,color='b',linestyle='--')
    ax.hist(np.log10(MV),mass_bins,histtype='step',density=b,color='g',label='vdXsec')

    if not b:
        ax.semilogy()
        fname = 'LogScale'
    else:
        fname = 'Normalized'
    ax.legend(loc='upper right',prop={'size':12})
    f.savefig(f'../Plots/MassHistogramComparison.{fname}.png',bbox_inches='tight',pad_inches=.1)