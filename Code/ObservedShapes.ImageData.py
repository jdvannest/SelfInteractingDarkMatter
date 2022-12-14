import pynbody,pickle,sys,pymp,warnings,argparse,os
import numpy as np
from pynbody.derived import rxy,az
import matplotlib.pylab as plt 
from math import pi,degrees
from scipy.optimize import curve_fit
from pynbody.plot.sph import image
from numpy.linalg import eig, inv
from matplotlib.patches import Ellipse
from numpy import sin,cos
plt.rcParams.update({'text.usetex':False})
#Sersic Function which is fit to the Surface-Brightness Pofiles of galaxy
def sersic(r, mueff, reff, n):
    return mueff + 2.5*(0.868*n-0.142)*((r/reff)**(1./n) - 1)
def myprint(string,clear=False):
    if clear:
        sys.stdout.write("\033[F")
        sys.stdout.write("\033[K") 
    print(string)
warnings.filterwarnings("ignore")

simpath = '/home/vannest/dwarf_volumes/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
HaloList = [1,2,3,5,6,7,10,11,13,14,24]
numproc = 10
#outfile = '../DataFiles/ObservedShapes.ImageData.pickle'
outfile = '../DataFiles/ObservedShapes.ImageData.pickle'
Data = {}

print(f'Loading sim...')
#Create shared dictionary for multiprocessing data and shared 1-element
# array to act as scalar for traking progress
SimData = pymp.shared.dict()
prog=pymp.shared.array((1,),dtype=int)
#Load in simulation
s = pynbody.load(simpath)
s.physical_units()
h = s.halos()
myprint(f'sim loaded.',clear=True)
print('\tWriting: 0%')
#Begin parallelized analysis
with pymp.Parallel(numproc) as pl:
    for i in pl.xrange(len(HaloList)):
        #Load in specific halo
        halonum = HaloList[i]
        halo = h[halonum]

        current={}

        #Imaging and Fitting Function
        def FitImage(centered_halo):
            #Create a smoothed SB profile to determine Reff and SB @ 2Reff
            try:
                Rvir = pynbody.analysis.halo.virial_radius(centered_halo)
            except:
                Rvir = 10
            prof = pynbody.analysis.profile.Profile(centered_halo.s,type='lin',min=.25,max=Rvir,ndim=2,nbins=int((Rvir-0.25)/0.1))
            #gband = gbandprofile(p['sb,b'],p['sb,v'])
            sb = prof['sb,v']
            #Smooth the sb profile
            smooth = np.nanmean(np.pad(sb.astype(float),(0,3-sb.size%3),mode='constant',constant_values=np.nan).reshape(-1,3),axis=1)
            try:
                y = smooth[:np.where(smooth>32)[0][0]+1]
            except:
                y = smooth
            #Create the x-data for the sb profile
            x = np.arange(len(y))*0.3 + 0.15
            x[0] = 0.05
            #Remove any NaN's from profile
            if True in np.isnan(y):
                x = np.delete(x,np.where(np.isnan(y)==True))
                y = np.delete(y,np.where(np.isnan(y)==True))
            #Perform Sersic fit on smoothed SB profile (if not all NaN's)
            if len(x)>0:
                #Initial guess values for curve_fit
                r0 = x[int(len(x)/2)]
                m0 = np.mean(y[:3])
                par,ign = curve_fit(sersic,x,y,p0=(m0,r0,1),bounds=([10,0,0.5],[40,100,16.5]))
                #find the array index closest to 2*Reff (Reff is par[1] from the Sersic function)
                ind = np.where(abs(prof['rbins']-2*par[1])==min(abs(prof['rbins']-2*par[1])))[0][0]
                #Get the v-band luminosity density at 2*Reff
                v = prof['v_lum_den'][ind]
                #Create a stellar particle filter around the halo and generate image
                width = 12*par[1]
                sphere = pynbody.filt.Sphere(width*np.sqrt(2)*1.01)
                im = image(s[sphere].s,qty='v_lum_den',width=width,units='kpc^-2',resolution=1000,show_cbar=False)
                return([v,im])            
            else:
                return([np.nan,np.nan])

        #Load in the halo, and use imaging/fitting function
        if len(halo.s)>0:
            pynbody.analysis.angmom.faceon(halo)
            faceSB,faceIm = FitImage(halo)
            pynbody.analysis.angmom.sideon(halo)
            sideSB,sideIm = FitImage(halo)
            current['SB_Faceon'] = faceSB
            current['Im_Faceon'] = faceIm
            current['SB_Sideon'] = sideSB
            current['Im_Sideon'] = sideIm
        else:
            current['SB_Faceon'] = np.nan
            current['Im_Faceon'] = np.nan
            current['SB_Sideon'] = np.nan
            current['Im_Sideon'] = np.nan

        #Print progress to terminal
        myprint(f'\tWriting: {round(float(prog[0]+1)/float(len(HaloList))*100,2)}%',clear=True)
        SimData[str(halonum)] = current
        prog[0]+=1
        del current

#Transfer the data from the shared dictionary to the main Data File (saving the shared
#dictionary just results in garbage data for some reason)
for halo in HaloList:
    Data[str(halo)] = SimData[str(halo)]

#Update the Data File
out = open(outfile,'wb')
pickle.dump(Data,out)
out.close()
print(f'File written.')