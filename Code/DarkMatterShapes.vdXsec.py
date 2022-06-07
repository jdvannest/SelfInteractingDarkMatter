import pickle,pymp,pynbody,sys,warnings
import numpy as np
from pynbody.analysis.halo import halo_shape
warnings.filterwarnings("ignore")
def myprint(string,clear=False):
    if clear:
        sys.stdout.write("\033[F")
        sys.stdout.write("\033[K") 
    print(string)
output_path = "../DataFiles/"
num_proc = 10

git = '/myhome2/users/vannest/SelfInteractingDarkMatter/'
simpath = '/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmbvdXsec.65536/storm.cosmo25cmbvdXsec.65536.065536'
filename = f'{git}DataFiles/DarkMatterShapes.Top200.vdXsec.pickle'

cont_frac = np.load(f'{git}DataFiles/ContaminationFraction.storm.vdXsec.z0.npy')
halo_list = np.where(cont_frac<.1)[0][:200]+1

print('Loading simulation...')
with open(f'{simpath}.z0.000.AHF_halos') as f:
    stat = f.readlines()
    del stat[0]
s = pynbody.load(simpath)
s.physical_units()
h = s.halos(dosort=True)
myprint('Simulation loaded.',clear=True)

print('Writing Pynbody Data: 0.00%')
SharedData = pymp.shared.dict()
prog=pymp.shared.array((1,),dtype=int)
with pymp.Parallel(num_proc) as pl:
    for i in pl.xrange(len(halo_list)):
        hid = halo_list[i]
        current = {}
        try:
            pynbody.analysis.angmom.faceon(h[hid])
            r,ba,ca,angle,Es = halo_shape(h[hid])
            current['b_pyn'] = ba
            current['c_pyn'] = ca
            current['rbins'] = r
        except:
            current['b_pyn'] = [np.nan]
            current['c_pyn'] = [np.nan]
            current['rbins'] = [np.nan]
        SharedData[str(hid)] = current
        prog[0]+=1
        myprint(f'Writing Pynbody Data: {round(prog[0]/len(halo_list)*100,2)}%',clear=True)

Data = {}
for halo in SharedData:
    Data[halo] = {}
    for key in ['b_pyn','c_pyn','rbins']:
        Data[halo][key] = SharedData[halo][key]

out = open(filename,'wb')
pickle.dump(Data,out)
out.close()
myprint('Pynbody Data Written.',clear=True)
