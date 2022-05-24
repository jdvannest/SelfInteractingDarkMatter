import argparse,sys
import numpy as np
def myprint(string,clear=False):
    if clear:
        sys.stdout.write("\033[F")
        sys.stdout.write("\033[K") 
    print(string)

parser = argparse.ArgumentParser(description='', usage='')
parser.add_argument('-c','--cross_section',choices=['CDM','SI3','SI10','vdXsec'],required=True)
parser.add_argument('-n','--npart',type=int,default=300)
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

print('Finding Halos')
Nhalo = 0
with open(AHFs[args.cross_section]) as f:
    AHF = f.readlines()
    del AHF[0]
index = 4 if args.cross_section=='vdXsec' else 3
for line in AHF:
    if int(line.split('\t')[index])>(args.npart-1): Nhalo+=1
myprint(f'{Nhalo} Halos Found.',clear=True)

Mvir,Rvir = np.zeros(Nhalo),np.zeros(Nhalo)
m_index = 3 if args.cross_section=='vdXsec' else 2
r_index = 11 if args.cross_section=='vdXsec' else 10
for i in np.arange(Nhalo):
    Mvir[i] = float(AHF[i].split('\t')[m_index])
    Rvir[i] = float(AHF[i].split('\t')[r_index])

np.save(f'../DataFiles/Mvir.N{args.npart}.{args.cross_section}.z0.npy',Mvir)
np.save(f'../DataFiles/Rvir.N{args.npart}.{args.cross_section}.z0.npy',Rvir)
print('Done')
