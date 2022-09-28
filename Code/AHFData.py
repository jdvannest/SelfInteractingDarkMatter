import argparse,os,pickle,sys
import numpy as np
def myprint(string,clear=False):
    if clear:
        sys.stdout.write("\033[F")
        sys.stdout.write("\033[K") 
    print(string)

parser = argparse.ArgumentParser(description='', usage='')
parser.add_argument('-s','--simulation',choices=['storm','h148'],required=True)
parser.add_argument('-c','--cross_section',choices=['CDM','SI3','SI10','SI50','vdXsec','VTS'],required=True)
parser.add_argument('-z','--redshift',choices=['z0','z0.04','z1','z2','z3','z4','z10','z13.5'],required=True)
#parser.add_argument('-n','--npart',type=int,default=64)
args = parser.parse_args()

#Update paths depending on machine
os.system(f'python Config.py')
config = pickle.load(open('Config.pickle','rb'))
AHF = config[args.simulation][args.redshift]['AHFs'][args.cross_section]

print('Finding Halos')
Nhalo = 0
with open(AHF) as f:
    AHF = f.readlines()
    del AHF[0]
#index = 4 
#for line in AHF:
#    if int(line.split('\t')[index])>(args.npart-1): Nhalo+=1
Nhalo = len(AHF)
myprint(f'{Nhalo} Halos Found.',clear=True)

Mvir,Rvir,cont,npart = np.zeros(Nhalo),np.zeros(Nhalo),np.zeros(Nhalo),np.zeros(Nhalo)
for i in np.arange(Nhalo):
    Mvir[i] = float(AHF[i].split('\t')[3])
    Rvir[i] = float(AHF[i].split('\t')[11])
    cont[i] = float(AHF[i].split('\t')[37])
    npart[i] = int(AHF[i].split('\t')[4])
cont = 1-cont

np.save(f'../DataFiles/Mvir.{args.simulation}.{args.cross_section}.{args.redshift}.npy',Mvir)
np.save(f'../DataFiles/Rvir.{args.simulation}.{args.cross_section}.{args.redshift}.npy',Rvir)
np.save(f'../DataFiles/Npart.{args.simulation}.{args.cross_section}.{args.redshift}.npy',npart)
np.save(f'../DataFiles/ContaminationFraction.{args.simulation}.{args.cross_section}.{args.redshift}.npy',cont)
print('Done')