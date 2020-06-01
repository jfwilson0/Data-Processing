# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 15:50:24 2017

@author: coke
"""

import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from decimal import *
getcontext().prec = 18
import matplotlib.colors as colors
import matplotlib.cm as cmx
from scipy.optimize import curve_fit
from collections import Counter
#plt.close('all')
c=0
a=10
j=199
out_str=''
name='ASD_00000'
ext='.dat'
n=np.zeros(j-a)
Rg=np.zeros(j-a)
Rgps=np.zeros(j-a)
Ai=np.zeros(j-a)
Bi=np.zeros(j-a)
t=np.zeros(j-a)
pk=np.zeros(j-a)

def func(x, a, b):
    return 10**(a*np.log(x)+b)
    #return a * np.exp(b * x)

def find10(m):
    mat=np.zeros((len(m), 3))
    mp = sorted(m, key=lambda a_entry: a_entry[1])
    
    for i in range(0,len(mp)):
        mat[i,0]=mp[i][0]
        mat[i,1]=mp[i][1]
        mat[i,2]=mp[i][2]

    nagg=len(mat)
    nu=int(round(nagg*0.1))
    o=np.zeros((nu, 3))
    mat=mat[::-1]
    
    for i in range(0,nu):
        o[i,:]=mat[i,:]
        
    return o
    
    
for i in range(a,j):
    #if (i%STEP==0):
        if i<1:
            out_str = name+ '0000'+ str(i) +ext            
        elif i<10:
            out_str = name+ '000'+ str(i) +ext
        elif i<100:
            out_str = name+ '00'+ str(i) +ext
        elif i<1000:
            out_str = name+ '0'+ str(i) +ext
        elif i<10000:
            out_str = name + str(i) +ext
        d = np.loadtxt(out_str, skiprows=0)
       
        matriz=find10(d)
        calculo=matriz[:,2]/matriz[:,1]
        pack=np.average(calculo)
        pk[c]=pack
#        ind = np.squeeze(np.asarray(d[:,1]))>1
#        j2=d[ind,:]
        Rgk=d[:,2]
        k=d[:,1]
        pr=1.0e-7
        one=np.ones(len(k))
        nps=3000-np.dot(k,one)
        Rgp=(0.775)*(pr)
        npp=np.ones(int(nps))
        Rgps=npp*Rgp
        kk=np.append(k,npp)
        Rgkk=np.append(Rgk,Rgps)
        N=len(k)
        one=np.ones(len(kk))
        A=(kk*kk)*(Rgkk*Rgkk)        
        B=kk*kk
        a=np.dot(A,one)
        b=np.dot(B,one)
        Rgpop=np.sqrt(a/b)
        t[c]=c
#        maximo=np.max(d[:,1])

#        #j2 = [i for i in d[:,1] if i >= (maximo/2)]
#        nm=np.average(j2[:,1])
#        Rgm=np.average(j2[:,2])        
        Rg[c]=Rgpop
        
        Ai[c], Bi[c] = np.polyfit(np.log(Rgk/pr), np.log(k), 1)
                    
        print(str(i) +' '+'first fit line:\ny = {:.2f} , {:.2f}'.format(Ai[c], Bi[c]))
        c=c+1
        
Rgtotal=np.average(Rg)        
Sigmd=np.std(Rg) 
        
plt.figure(1,facecolor="white")
ax2 = plt.subplot(111)
ax2.plot(t,Rg,'o',c='m',label='',linewidth=2)
ax3 = ax2.twinx()
ax3.plot(t,Ai,'o',c='r',label='',linewidth=2)
ax3.set_ylabel('df', color='k',fontsize=20)
ax2.set_xlabel('Time (-)',fontsize=20)
ax2.set_ylabel('Rg', color='k',fontsize=20)
plt.legend(loc='lower center', ncol=3,fancybox=True, shadow=True, fontsize=12)
ax2.tick_params(axis='y', labelsize=16)
ax2.tick_params(axis='x', labelsize=16)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_yscale('log')
plt.grid()
plt.show()


RG=Rgk/pr

Rg_avg = np.average(Rg[120:190])
print(Rg_avg)

np.savetxt('Rgdf16004DMTSEvE1.txt',np.transpose([t,Rg,Ai,pk]),delimiter=' ')


####################################################################################



popt2=np.polyfit(np.log(RG), np.log(k),1)
p = np.poly1d(popt2)

erreg=[0,1,2,3,4,5]


plt.figure(2,facecolor="white")
ax2 = plt.subplot(111)
plt.plot(np.log(RG), np.log(k),'o')
ax2.plot(erreg,p(erreg),c='r') # polyfit del ultimo archivo
plt.grid()
plt.show()


plt.figure(3,facecolor="white")
ax2 = plt.subplot(111)
ax2.plot(t,1/pk,'o',c='b',label='',linewidth=2)
ax2.set_xlabel('Time (-)',fontsize=20)
ax2.set_ylabel('Packing number', color='k',fontsize=20)
plt.legend(loc='lower center', ncol=3,fancybox=True, shadow=True, fontsize=12)
ax2.tick_params(axis='y', labelsize=16)
ax2.tick_params(axis='x', labelsize=16)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.grid()
plt.show()