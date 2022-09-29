import pickle
import numpy as np
from math import pi

s = pickle.load(open('../DataFiles/DarkMatterShapes.Top200.vdXsec.pickle','rb'))
eps = []

def vol(a,b,c):
    return((4/3)*pi*a*b*c)

for h in s:
    a,ba,ca,n = s[h]['rbins'],s[h]['b_pyn'],s[h]['c_pyn'],s[h]['n']
    i= np.argmin(np.array(s[h]['rbins'] - (s[h]['rbins'][-1]*.1 )))
    try:
        eps.append( (vol(a[i],a[i]*ba[i],a[i]*ca[i])-vol(a[i-1],a[i-1]*ba[i-1],a[i-1]*ca[i-1]))/np.sqrt(n[i]) * np.abs((n[i-1]/(vol(a[i-1],a[i-1]*ba[i-1],a[i-1]*ca[i-1])-vol(a[i-2],a[i-2]*ba[i-2],a[i-2]*ca[i-2])))-(n[i+1]/(vol(a[i+1],a[i+1]*ba[i+1],a[i+1]*ca[i+1])-vol(a[i],a[i]*ba[i],a[i]*ca[i])))) )
    except:
        eps.append(np.NaN)
print(eps)