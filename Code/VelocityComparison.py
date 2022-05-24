import numpy as np
import matplotlib.pylab as plt

for z in ['.z0','.z1']:
    M0 = np.load(f'../DataFiles/Mvir.CDM{z}.npy')
    R0 = np.load(f'../DataFiles/Rvir.CDM{z}.npy')
    M3 = np.load(f'../DataFiles/Mvir.SI3{z}.npy')
    R3 = np.load(f'../DataFiles/Rvir.SI3{z}.npy')
    M10 = np.load(f'../DataFiles/Mvir.SI10{z}.npy')
    R10 = np.load(f'../DataFiles/Rvir.SI10{z}.npy')
    MV = np.load(f'../DataFiles/Mvir.vdXsec{z}.npy')
    RV = np.load(f'../DataFiles/Rvir.vdXsec{z}.npy')

    v_bins = np.logspace(0,2.65,100)
    G = 4.30091e-3 #pc Msol^-1 (km/s)^2

    f,ax = plt.subplots(1,1,figsize=(6,4))
    ax.set_xlabel(r'V$_{vir}$ [km/s]',fontsize=15)
    ax.set_ylabel('N',fontsize=15)
    ax.semilogy()
    ax.semilogx()

    ax.hist(np.sqrt((M0*G)/(R0*1000)),v_bins,histtype='step',color='k',label='CDM')
    ax.hist(np.sqrt((M3*G)/(R3*1000)),v_bins,histtype='step',color='r',label='SI3')
    ax.hist(np.sqrt((M10*G)/(R10*1000)),v_bins,histtype='step',color='b',label='SI10')
    ax.hist(np.sqrt((MV*G)/(RV*1000)),v_bins,histtype='step',color='g',label='vdXsec')

    ax.legend(loc='upper right',prop={'size':12})
    f.savefig(f'../Plots/VelocityComparison{z}.png',bbox_inches='tight',pad_inches=.1)

    #Normalized
    f,ax = plt.subplots(1,1,figsize=(6,4))
    ax.set_xlabel(r'V$_{vir}$ [km/s]',fontsize=15)
    ax.set_ylabel('Normalized',fontsize=15)
    ax.semilogx()

    ax.hist(np.sqrt((M0*G)/(R0*1000)),v_bins,histtype='step',density=True,color='k',label='CDM')
    ax.hist(np.sqrt((M3*G)/(R3*1000)),v_bins,histtype='step',density=True,color='r',label='SI3')
    ax.hist(np.sqrt((M10*G)/(R10*1000)),v_bins,histtype='step',density=True,color='b',label='SI10')
    ax.hist(np.sqrt((MV*G)/(RV*1000)),v_bins,histtype='step',density=True,color='g',label='vdXsec')

    ax.legend(loc='upper right',prop={'size':12})
    f.savefig(f'../Plots/VelocityComparison{z}.Normalized.png',bbox_inches='tight',pad_inches=.1)