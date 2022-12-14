import pickle
import numpy as np
from math import pi

s = pickle.load(open('../DataFiles/DarkMatterShapes.Top200.vdXsec.pickle','rb'))
eps = []

def vol(a,b,c):
    return((4/3)*pi*a*b*c)

for h in s:
    a,ba,ca,n = s[h]['rbins'],s[h]['b_pyn'],s[h]['c_pyn'],s[h]['n']
    b,c = a*ba,a*ca
    v = np.zeros(len(a))
    for i in np.arange(len(a)):
        if i==0: v[i] = vol(a[i],b[i],c[i])
        else: v[i] = vol(a[i],b[i],c[i])-vol(a[i-1],b[i-1],c[i-1])
    i= np.argmin(abs(s[h]['rbins'] - (s[h]['rbins'][-1]*.1 )))
    eps.append( v[i]/np.sqrt(n[i]) * np.abs( n[i-1]/v[i-1] - n[i+1]/v[i+1] ) )


h = '200'
a,ba,ca,n = s[h]['rbins'],s[h]['b_pyn'],s[h]['c_pyn'],s[h]['n']
b,c = a*ba,a*ca
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pylab as plt
f = plt.figure()
ax = f.add_subplot(111, projection='3d')

if True:
    for i in np.arange(len(a)):
        #coefs = (a[i],b[i],c[i])
        #rx, ry, rz = 1/np.sqrt(coefs)
        u = np.linspace(0,2*pi,100)
        v = np.linspace(0,pi,100)
        x = a[i] * np.outer(np.cos(u),np.sin(v))
        y = b[i] * np.outer(np.sin(u),np.sin(v))
        z = c[i] * np.outer(np.ones_like(u),np.cos(v))
        ax.plot_surface(x,y,z,rstride=4,cstride=4,color='b',alpha=.3)
        #for axis in 'xyz': getattr(ax, 'set_{}lim'.format(axis))((-a[-1],a[-1]))
        ax.set_xlim([-a[-1],a[-1]])
        ax.set_ylim([-a[-1],a[-1]])
        ax.set_zlim([-a[-1],a[-1]])
        f.savefig(f'/Users/jdvannest/Desktop/Ellipsoid/Ell.{i:03d}.png',bbox_inches='tight',pad_inches=.1)
        ax.clear()

import os,imageio
imagenames = os.listdir('/Users/jdvannest/Desktop/Ellipsoid/')
imagenames = np.sort(np.array(imagenames))
images = []
for name in imagenames:
    images.append(imageio.imread(f'/Users/jdvannest/Desktop/Ellipsoid/{name}'))
imageio.mimsave('/Users/jdvannest/Desktop/Ellipsoid/AAA.gif', images, duration=.1)