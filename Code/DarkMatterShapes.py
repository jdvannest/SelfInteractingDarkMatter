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
parser.add_argument("-c","--cross_section",required=True,choices=['3','10','30','50'])
parser.add_argument("-p","--pynbody",action="store_true")
args = parser.parse_args()

Data = {}

print('Loading simulation...')
simpath = f'/data2/akaxia/storm.cosmo25cmbSI{args.cross_section}.4096/storm.cosmo25cmbSI{args.cross_section}.4096.004096'
with open(f'{simpath}.0000.z0.000.AHF_halos') as f:
    stat = f.readlines()
s = pynbody.load(simpath)
s.physical_units()
h = s.halos(dosort=True)
myprint('Simulation loaded.\nWriting AHF Data',clear=True)

#stat npart:4 , b:24 , c:25
resolved,hid,halos = [True,1,[]]
while resolved:
    if int(stat[hid].split()[4]) < 1000:
        resolved = False
    else:
        Data[str(hid)] = {'b':np.nan,'c':np.nan,'b_pyn':[np.nan],'c_pyn':[np.nan],'rbins':[np.nan]}
        Data[str(hid)]['b'] = float(stat[hid].split()[24])
        Data[str(hid)]['c'] = float(stat[hid].split()[25])
        halos.append(hid)
        hid+=1

out = open(f'{output_path}DarkMatterShapes.SI{args.cross_section}.pickle','wb')
pickle.dump(Data,out)
out.close()
myprint('AHF Data Written.',clear=True)

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
            myprint(f'Writing Pynbody Data: {round(hid/len(halos),2)}%',clear=True)
            prog+=1
    
    for halo in Data:
        for key in ['b_pyn','c_pyn','rbins']:
            Data[halo][key] = SharedData[halo][key]
    
    out = open(f'{output_path}DarkMatterShapes.SI{args.cross_section}.pickle','wb')
    pickle.dump(Data,out)
    out.close()
    myprint('Pynbody Data Written.',clear=True)