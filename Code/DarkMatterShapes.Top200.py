import argparse,pickle,pymp,pynbody,sys,warnings
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


parser = argparse.ArgumentParser(description="Calculates Shapes of dark matter halos"
                                +"using pynbody's built in shape function", 
                                usage="DarkMatterShapes.py -c 3")
parser.add_argument("-c","--cross_section",required=True,choices=['cdm','3','10','30','50'])
parser.add_argument("-p","--pynbody",action="store_true")
args = parser.parse_args()

Data = {}
if args.cross_section == 'cdm':
    simpath = '/data/REPOSITORY/dwarf_volumes/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.004096'
    filename = f'{output_path}DarkMatterShapes.Top200.CDM.pickle'
else:
    simpath = f'/data2/akaxia/storm.cosmo25cmbSI{args.cross_section}.4096/storm.cosmo25cmbSI{args.cross_section}.4096.004096'
    filename = f'{output_path}DarkMatterShapes.Top200.SI{args.cross_section}.pickle'

with open('/data2/akaxia/storm_uncontam_halos_top200.txt') as f:
    L = f.readlines()
halo_list = [int(float(i.rstrip('\n'))) for i in L]

print('Loading simulation...')
with open(f'{simpath}.0000.z0.000.AHF_halos') as f:
    stat = f.readlines()
s = pynbody.load(simpath)
s.physical_units()
h = s.halos(dosort=True)
myprint('Simulation loaded.\nWriting AHF Data',clear=True)

#stat npart:4 , b:24 , c:25
for hid in halo_list:
    Data[str(hid)] = {'b':np.nan,'c':np.nan,'b_pyn':[np.nan],'c_pyn':[np.nan],'rbins':[np.nan]}
    Data[str(hid)]['b'] = float(stat[hid].split()[24])
    Data[str(hid)]['c'] = float(stat[hid].split()[25])


out = open(filename,'wb')
pickle.dump(Data,out)
out.close()
myprint(f'AHF Data Written. {len(halo_list)} Resolved Halos Found.',clear=True)

if args.pynbody:
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
    
    for halo in SharedData:
        for key in ['b_pyn','c_pyn','rbins']:
            Data[halo][key] = SharedData[halo][key]
    
    out = open(filename,'wb')
    pickle.dump(Data,out)
    out.close()
    myprint('Pynbody Data Written.',clear=True)