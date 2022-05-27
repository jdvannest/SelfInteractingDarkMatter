import argparse
import numpy as np
import matplotlib.pylab as plt

parser = argparse.ArgumentParser()
parser.add_argument('-z','--redshift',choices=['z0','z1','z2','z3','z4'],required=True)
parser.add_argument('-n','--npart',default=64)
args = parser.parse_args()

#Load Data
M0 = np.load(f'../DataFiles/Mvir.N{args.npart}.CDM.{args.redshift}.npy')
R0 = np.load(f'../DataFiles/Rvir.N{args.npart}.CDM.{args.redshift}.npy')
M3 = np.load(f'../DataFiles/Mvir.N{args.npart}.SI3.{args.redshift}.npy')
R3 = np.load(f'../DataFiles/Rvir.N{args.npart}.SI3.{args.redshift}.npy')
M10 = np.load(f'../DataFiles/Mvir.N{args.npart}.SI10.{args.redshift}.npy')
R10 = np.load(f'../DataFiles/Rvir.N{args.npart}.SI10.{args.redshift}.npy')
MV = np.load(f'../DataFiles/Mvir.N{args.npart}.vdXsec.{args.redshift}.npy')
RV = np.load(f'../DataFiles/Rvir.N{args.npart}.vdXsec.{args.redshift}.npy')
#Apply contamination fraction limit
M0c = np.load(f'../DataFiles/ContaminationFraction.N{args.npart}.CDM.{args.redshift}.npy')
M3c = np.load(f'../DataFiles/ContaminationFraction.N{args.npart}.SI3.{args.redshift}.npy')
M10c = np.load(f'../DataFiles/ContaminationFraction.N{args.npart}.SI10.{args.redshift}.npy')
MVc = np.load(f'../DataFiles/ContaminationFraction.N{args.npart}.vdXsec.{args.redshift}.npy')
contam_limit=.1
M0 = M0[M0c<contam_limit]
R0 = R0[M0c<contam_limit]
M3 = M3[M3c<contam_limit]
R3 = R3[M3c<contam_limit]
M10 = M10[M10c<contam_limit]
R10 = R10[M10c<contam_limit]
MV = MV[MVc<contam_limit]
RV = RV[MVc<contam_limit]

v_bins = np.logspace(-.3,2.65,100)
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
f.savefig(f'../Plots/VelocityComparison.N{args.npart}.{args.redshift}.png',bbox_inches='tight',pad_inches=.1)

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
f.savefig(f'../Plots/VelocityComparison.N{args.npart}.{args.redshift}.Normalized.png',bbox_inches='tight',pad_inches=.1)