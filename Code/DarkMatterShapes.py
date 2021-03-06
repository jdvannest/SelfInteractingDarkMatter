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
    filename = f'{output_path}DarkMatterShapes.CDM.pickle'
else:
    simpath = f'/data2/akaxia/storm.cosmo25cmbSI{args.cross_section}.4096/storm.cosmo25cmbSI{args.cross_section}.4096.004096'
    filename = f'{output_path}DarkMatterShapes.SI{args.cross_section}.pickle'

print('Loading simulation...')
with open(f'{simpath}.0000.z0.000.AHF_halos') as f:
    stat = f.readlines()
s = pynbody.load(simpath)
s.physical_units()
h = s.halos(dosort=True)
myprint('Simulation loaded.\nWriting AHF Data',clear=True)

#stat npart:4 , b:24 , c:25
Data['1'] = {'b':float(stat[1].split()[24]),'c':float(stat[1].split()[25]),'b_pyn':[np.nan],'c_pyn':[np.nan],'rbins':[np.nan]}
resolved_num,hid,halos,min_dm_mass = [True,2,[1],8067.581622282078]
while resolved_num:
    if int(stat[hid].split()[4]) < 1000:
        resolved_num = False
    else:
        if len(h[hid].dm['mass'][h[hid].dm['mass']>min_dm_mass])/len(h[hid].dm['mass']) < 0.1:
            Data[str(hid)] = {'b':np.nan,'c':np.nan,'b_pyn':[np.nan],'c_pyn':[np.nan],'rbins':[np.nan]}
            Data[str(hid)]['b'] = float(stat[hid].split()[24])
            Data[str(hid)]['c'] = float(stat[hid].split()[25])
            halos.append(hid)
        hid+=1

out = open(filename,'wb')
pickle.dump(Data,out)
out.close()
myprint(f'AHF Data Written. {len(halos)} Resolved Halos Found.',clear=True)

if args.pynbody:
    print('Writing Pynbody Data: 0.00%')
    SharedData = pymp.shared.dict()
    prog=pymp.shared.array((1,),dtype=int)
    with pymp.Parallel(num_proc) as pl:
        for i in pl.xrange(len(halos)):
            hid = halos[i]
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
            myprint(f'Writing Pynbody Data: {round(prog[0]/len(halos)*100,2)}%',clear=True)
    
    for halo in Data:
        for key in ['b_pyn','c_pyn','rbins']:
            Data[halo][key] = SharedData[halo][key]
    
    out = open(filename,'wb')
    pickle.dump(Data,out)
    out.close()
    myprint('Pynbody Data Written.',clear=True)