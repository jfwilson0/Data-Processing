# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 12:09:56 2014

@author: coke
"""
import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from decimal import *
getcontext().prec = 18
import matplotlib.colors as colors
import matplotlib.cm as cmx


######################################################################################################
def overlap(d,f):
    Npart=3000       #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    h_max=7.3204406453619177E-010/2    #!!!!!!!!!!!!!!!!!!!!!!
    psize=1.0e-7
    Rlength=2*psize+2*h_max
    ID1=d[:,0];ID2=d[:,1];x1=d[:,2];y1=d[:,3];z1=d[:,4];x2=d[:,5];y2=d[:,6];z2=d[:,7]
    d_x=x2-x1;d_y=y2-y1;d_z=z2-z1
    FN=np.sqrt(d[:,8]*d[:,8]+d[:,9]*d[:,9]+d[:,10]*d[:,10])
    FT=np.sqrt(d[:,11]*d[:,11]+d[:,12]*d[:,12]+d[:,13]*d[:,13])
    l=np.zeros(len(x1))
    contador=0
    
    
    for i in range(0, len(x1)):
        l[i]=np.sqrt(d_x[i]*d_x[i]+d_y[i]*d_y[i]+d_z[i]*d_z[i])
        if l[i]<(Rlength):
            contador=contador+1
    O=0
    contactno=len(x1)
    
    
    overeal=np.zeros(contador)
    FrealN=np.zeros(contador)
    FrealT=np.zeros(contador)
    co=0 
    for i in range(0, len(d_x)):
        if l[i]<(Rlength):
            overeal[co]=l[i]
            FrealN[co]=FN[i]
            FrealT[co]=FT[i]
            co=co+1

    Aov=np.mean(overeal)
    O=(Rlength-Aov)*1e9
    AFN=np.mean(FrealN)*1000
    AFT=np.mean(FrealT)*1000
    return O,contactno, AFT, AFN
######################################################################################################

STEP=1
STOPT=199
SR=108         #!!!!!!!!!!!!!!
printIt=5000   #!!!!!!!!!!!!!!
K=0.4            #!!!!!!!!!!!!!!
delta='5.0'          #!!!!!!!!!!!!!!
timestep=1e-10

Nm=0
for i in range(STEP,STOPT):
    if (i%STEP==0):
        Nm=Nm+1
cuenta=0
out_str=''
out_str2=''
name='contacts_00000'
name2='elems_00000'
ext='.dat'
nueve='9'
S1=np.zeros(Nm+1);S2=np.zeros(Nm+1);S3=np.zeros(Nm+1);S4=np.zeros(Nm+1);S5=np.zeros(Nm+1)
I=np.zeros(Nm+1)
for i in range(STEP,STOPT):
    if (i%STEP==0):
        if i<1:
            out_str = name+ '0000'+ `i` +ext
            out_str2 = name2+ '0000'+ `i` + ext  
        elif i<10:
            out_str = name+ '000'+ `i` +ext
            out_str2 = name2+ '000'+ `i` + ext
        elif i<100:
            out_str = name+ '00'+`i` +ext
            out_str2 = name2+ '00'+`i` + ext
        elif i<1000:
            out_str = name+ '0'+`i` +ext
            out_str2 = name2+ '0'+`i` + ext
        elif i<10000:
            out_str = name +`i` +ext
            out_str2 = name2 +`i` + ext
        d = np.loadtxt(out_str, skiprows=0)
        f = np.loadtxt(out_str2, skiprows=0)
        cuenta=cuenta+1
        S1[cuenta],_,_,_=overlap(d,f)
        _,S2[cuenta],_,_=overlap(d,f)
        _,_,S3[cuenta],_=overlap(d,f)
        _,_,_,S4[cuenta]=overlap(d,f)
        I[cuenta]=i*printIt*SR*1e4*timestep
    


plt.figure(1,facecolor="white")
ax2 = plt.subplot(111)
ax2.plot(I,S4,':o',c='m',label='Average Normal force')
ax2.plot(I,S3,':*',c='m',label='Average Tangential force')
ax2.set_xlabel('Time',fontsize=16)
ax2.set_ylabel('Force [N]', color='k',fontsize=16)
plt.legend(scatterpoints=1, loc='upper left', ncol=1, fontsize=12)
ax4= plt.twinx()
p2=ax4.plot(I,S1,'-s',c='k',label='Average overlap')
plt.legend(scatterpoints=1, loc='upper right', ncol=1, fontsize=12)
#ax2.set_ylim(0, 0.000014)
#ax4.set_ylim(0,1.4)
plt.show()


plt.figure(2,facecolor="white")
ax2 = plt.subplot(111)
ax2.plot(I,S2,'--',c='b',label='Number of contacts E05')
plt.legend(scatterpoints=1, loc='upper right', ncol=1, fontsize=12)
ax3 = ax2.twinx()
ax2.set_xlabel('Time',fontsize=16)
ax2.set_ylabel('Contact no', color='k',fontsize=16)
plt.legend(scatterpoints=1, loc='upper left', ncol=1, fontsize=12)
#ax2.set_ylim(0, 0.000014)
#ax4.set_ylim(0,1.4)
plt.show()

Overlap_ave = np.average(S1[120:190])
print(Overlap_ave)


np.savetxt('stats'+'_'+`SR`+'_'+delta+'_'+`K`+'DMTSEvE1.txt',np.transpose([I,S1,S2,S3,S4]),delimiter=' ')












