import argparse,pickle,pynbody,sys
import numpy as np
from pynbody.analysis.halo import halo_shape
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
args = parser.parse_args()

Data = {}

print('Loading simulation...')
simpath = f'/data2/akaxia/storm.cosmo25cmbSI{args.cross_section}.4096/storm.cosmo25cmbSI{args.cross_section}.4096.004096'
with open(f'{simpath}.0000.z0.000.AHF_halos') as f:
    stat = f.readlines()
s = pynbody.load(simpath)
s.physical_units()
h = s.halos(dosort=True)
myprint('Simulation loaded.\n')

#stat npart:4 , b:24 , c:25
resolved,hid = [True,1]
while resolved:
    myprint(f'Current DM Particle Count: {stat[hid].split()[4]}',clear=True)
    if int(stat[hid].split()[4]) < 1000:
        resolved = False
    else:
        Data[str(hid)] = {'b':np.nan,'c':np.nan,'b_pyn':[],'c_pyn':[],'rbins':[]}
        Data[str(hid)]['b'] = float(stat[hid].split()[24])
        Data[str(hid)]['c'] = float(stat[hid].split()[25])
        try:
            r,ba,ca,angle,Es = halo_shape(h[hid])
            Data[str(hid)]['b_pyn'] = ba
            Data[str(hid)]['c_pyn'] = ca
            Data[str(hid)]['rbins'] = r
        except:
            pass
        hid+=1

out = open(f'{output_path}DarkMatterShapes.SI{args.cross_section}.pickle','wb')
pickle.dump(Data,out)
out.close()
print('File Updated.')