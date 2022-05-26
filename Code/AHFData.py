import argparse,os,pickle,sys
import numpy as np
def myprint(string,clear=False):
    if clear:
        sys.stdout.write("\033[F")
        sys.stdout.write("\033[K") 
    print(string)

parser = argparse.ArgumentParser(description='', usage='')
parser.add_argument('-c','--cross_section',choices=['CDM','SI3','SI10','vdXsec'],required=True)
parser.add_argument('-z','--redshift',choices=['z0','z1','z2','z3','z4'],required=True)
parser.add_argument('-n','--npart',type=int,default=64)
args = parser.parse_args()

#Update paths depending on machine
os.system(f'python Config.py')
config = pickle.load(open('Config.pickle','rb'))
AHF = config[args.redshift]['AHFs'][args.cross_section]

print('Finding Halos')
Nhalo = 0
with open(AHF) as f:
    AHF = f.readlines()
    del AHF[0]
index = 4 
for line in AHF:
    if int(line.split('\t')[index])>(args.npart-1): Nhalo+=1
myprint(f'{Nhalo} Halos Found.',clear=True)

Mvir,Rvir,cont = np.zeros(Nhalo),np.zeros(Nhalo),np.zeros(Nhalo)
m_index = 3
r_index = 11
c_index = 37
for i in np.arange(Nhalo):
    Mvir[i] = float(AHF[i].split('\t')[m_index])
    Rvir[i] = float(AHF[i].split('\t')[r_index])
    cont[i] = float(AHF[i].split('\t')[c_index])

np.save(f'../DataFiles/Mvir.N{args.npart}.{args.cross_section}.{args.redshift}.npy',Mvir)
np.save(f'../DataFiles/Rvir.N{args.npart}.{args.cross_section}.{args.redshift}.npy',Rvir)
np.save(f'../DataFiles/ContaminationFraction.N{args.npart}.{args.cross_section}.{args.redshift}.npy',cont)
print('Done')
