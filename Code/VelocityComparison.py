import argparse
import numpy as np
import matplotlib.pylab as plt

parser = argparse.ArgumentParser()
parser.add_argument('-z','--redshift',choices=['z0','z1','z2','z3','z4'],required=True)
parser.add_argument('-n','--npart',default=64)
args = parser.parse_args()

#Load Data
N0 = np.load(f'../DataFiles/Npart.CDM.{args.redshift}.npy')
N3 = np.load(f'../DataFiles/Npart.SI3.{args.redshift}.npy')
N10 = np.load(f'../DataFiles/Npart.SI10.{args.redshift}.npy')
NV = np.load(f'../DataFiles/Npart.vdXsec.{args.redshift}.npy')
Nh = np.load(f'../DataFiles/Npart.h148.{args.redshift}.npy')
M0 = np.load(f'../DataFiles/Mvir.CDM.{args.redshift}.npy')
M3 = np.load(f'../DataFiles/Mvir.SI3.{args.redshift}.npy')
M10 = np.load(f'../DataFiles/Mvir.SI10.{args.redshift}.npy')
MV = np.load(f'../DataFiles/Mvir.vdXsec.{args.redshift}.npy')
Mh = np.load(f'../DataFiles/Mvir.h148.{args.redshift}.npy')
R0 = np.load(f'../DataFiles/Rvir.CDM.{args.redshift}.npy')
R3 = np.load(f'../DataFiles/Rvir.SI3.{args.redshift}.npy')
R10 = np.load(f'../DataFiles/Rvir.SI10.{args.redshift}.npy')
RV = np.load(f'../DataFiles/Rvir.vdXsec.{args.redshift}.npy')
Rh = np.load(f'../DataFiles/Rvir.h148.{args.redshift}.npy')
C0 = np.load(f'../DataFiles/ContaminationFraction.CDM.{args.redshift}.npy')
C3 = np.load(f'../DataFiles/ContaminationFraction.SI3.{args.redshift}.npy')
C10 = np.load(f'../DataFiles/ContaminationFraction.SI10.{args.redshift}.npy')
CV = np.load(f'../DataFiles/ContaminationFraction.vdXsec.{args.redshift}.npy')
Ch = np.load(f'../DataFiles/ContaminationFraction.h148.{args.redshift}.npy')
#Apply particle count limit
M0 = M0[N0>(int(args.npart)-1)]
C0 = C0[N0>(int(args.npart)-1)]
R0 = R0[N0>(int(args.npart)-1)]
M3 = M3[N3>(int(args.npart)-1)]
C3 = C3[N3>(int(args.npart)-1)]
R3 = R3[N3>(int(args.npart)-1)]
M10 = M10[N10>(int(args.npart)-1)]
C10 = C10[N10>(int(args.npart)-1)]
R10 = R10[N10>(int(args.npart)-1)]
MV = MV[NV>(int(args.npart)-1)]
CV = CV[NV>(int(args.npart)-1)]
RV = RV[NV>(int(args.npart)-1)]
Mh = Mh[Nh>(int(args.npart)-1)]
Ch = Ch[Nh>(int(args.npart)-1)]
Rh = Rh[Nh>(int(args.npart)-1)]
#Apply contamination fraction limit
contam_limit=.1
M0 = M0[C0<contam_limit]
M3 = M3[C3<contam_limit]
M10 = M10[C10<contam_limit]
MV = MV[CV<contam_limit]
Mh = Mh[Ch<contam_limit]
R0 = R0[C0<contam_limit]
R3 = R3[C3<contam_limit]
R10 = R10[C10<contam_limit]
RV = RV[CV<contam_limit]
Rh = Rh[Ch<contam_limit]

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
ax.hist(np.sqrt((Mh*G)/(Rh*1000)),v_bins,histtype='step',color='orange',label='h148')

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
ax.hist(np.sqrt((Mh*G)/(Rh*1000)),v_bins,histtype='step',density=True,color='orange',label='h148')

ax.legend(loc='upper right',prop={'size':12})
f.savefig(f'../Plots/VelocityComparison.N{args.npart}.{args.redshift}.Normalized.png',bbox_inches='tight',pad_inches=.1)