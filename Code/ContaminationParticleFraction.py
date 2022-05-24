import argparse,pynbody,sys
import numpy as np
def myprint(string,clear=False):
    if clear:
        sys.stdout.write("\033[F")
        sys.stdout.write("\033[K") 
    print(string)

parser = argparse.ArgumentParser(description='', usage='')
parser.add_argument('-c','--cross_section',choices=['CDM','SI3','SI10','vdXsec'],required=True)
parser.add_argument('-n','--npart',default=300)
args = parser.parse_args()

simpaths = {
    'CDM':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.004096',
    'SI3':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmbSI3.4096/storm.cosmo25cmbSI3.4096.004096',
    'SI10':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmbSI10.4096/storm.cosmo25cmbSI10.4096.004096',
    'vdXsec':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmbvdXsec.65536/storm.cosmo25cmbvdXsec.65536.065536'
}
AHFs = {
    'CDM':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.004096.0000.z0.000.AHF_halos',
    'SI3':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmbSI3.4096/storm.cosmo25cmbSI3.4096.004096.0000.z0.000.AHF_halos',
    'SI10':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmbSI10.4096/storm.cosmo25cmbSI10.4096.004096.0000.z0.000.AHF_halos',
    'vdXsec':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmbvdXsec.65536/storm.cosmo25cmbvdXsec.65536.065536.z0.000.AHF_halos'
}

print('Loading Sim...')
s = pynbody.load(simpaths[args.cross_section])
s.physical_units()
h = s.halos()
myprint('Sim Loaded.',clear=True)

print('Finding Halos')
Nhalo = 0
with open(AHFs[args.cross_section]) as f:
    AHF = f.readlines()
    del AHF[0]
index = 4 if args.cross_section=='vdXsec' else 3
for line in AHF:
    if int(line.split('\t')[index])>(args.npart-1): Nhalo+=1
myprint(f'{Nhalo} Halos Found.',clear=True)

prog,contam_frac,min_mass = 0,np.zeros(Nhalo),min(h[1].d['mass'])
print('Writing: 0.00%')
for i in np.arange(Nhalo):
    #Contamination Fraction is fraction of dm particles that are more massive than the
    #minimum dm particles mass (fraction of 'low res' particles in halo)
    contam_frac[i] = len(h[i+1].d['mass'][h[i+1].d['mass']>min_mass])/len(h[i+1].d['mass'])
    myprint(f'Writing: {round((i+1)/Nhalo*100,2)}%',clear=True)

np.save(f'../DataFiles/ContaminationParticleFraction.N{args.npart}.{args.cross_section}.z0.npy',contam_frac)
myprint('Done',clear=True)
