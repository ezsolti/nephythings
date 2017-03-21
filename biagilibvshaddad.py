"""this script performs the verificiation of Bolsig against Haddad Drift Velocity of 
Electrons in Nitrogen-Argon Mixtures

It requires bolsigminus.exe and a reference input file for that code.

Zs. Elter, 2017."""
import os
import math

def InputMaker(GasTemperature,ReducedElectricField,ArRat,NRat):
    bolsigRefFile = open('input-reference.dat','r')
    bolsigStr = bolsigRefFile.read()
    bolsigRefFile.close()

    bolsigStr = bolsigStr.replace("ReducedElectricField",str(ReducedElectricField))
    bolsigStr = bolsigStr.replace("GasTemp",str(GasTemperature))
    bolsigStr = bolsigStr.replace("ArRat",str(ArRat))
    bolsigStr = bolsigStr.replace("NRat",str(NRat))

    bolsigFinFile = open('bolsig-pyFC.dat','w')
    bolsigFinFile.write(bolsigStr)
    bolsigFinFile.close()
    
def Runner():
    os.system('bolsigminus.exe bolsig-pyFC.dat')
 
def OutputExtractor():
    bolsigOutFile = open('bolsigcalc.dat','r')

    mobilityFlag=False
    for line in bolsigOutFile:
        x=line.strip().split()
        if mobilityFlag:
            mobilityFlag=False                    
            electronRawMobility=float(x[1])
        if len(x)>=3 and x[2]=='Mobility': 
            mobilityFlag=True

    bolsigOutFile.close()
    return electronRawMobility

######################################################################
#Input
GasTemperature=293   #K

ArRat=[0.999,0.990,0.950]
NRat=[0.001,0.010,0.050]

#for i in range(len(ArRat)):
#    SumOc=ArRat[i]+2*NRat[i]
#    ArRat[i]=ArRat[i]/SumOc
#    NRat[i]=(NRat[i]*2)/SumOc

#reduced electric field in Haddad's paper
REr01=[0.1,0.12,0.14,0.17,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8,1.0,1.2,1.4,1.7,2.0,2.5,3.0,3.5,4.0,5.0]
REr1=[0.06,0.07,0.08,0.1,0.12,0.14,0.17,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8,1.0,1.2,1.4,1.7,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,10.0]
REr5=[0.1,0.12,0.14,0.17,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8,1.0,1.2,1.4,1.7,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,10.0]
#electron drift velocity in Haddad's paper
haddad01=[1.75,1.88,2.01,2.22,2.44,2.82,3.22,3.63,4.03,4.81,5.55,6.22,6.82,7.77,8.30,8.48,8.28,7.84,7.11,6.63,6.38,6.28,6.46]
haddad1=[2.08,2.12,2.15,2.24,2.35,2.48,2.69,2.90,3.28,3.68,4.08,4.49,5.31,6.14,6.98,7.81,9.48,11.1,12.7,14.9,17.0,19.8,21.6,22.4,22.5,21.3,20.0,19.1,18.6,18.3]
haddad5=[3.3,3.34,3.39,3.5,3.63,3.9,4.2,4.53,4.87,5.55,6.25,6.95,7.65,9.04,10.5,11.8,13.9,15.9,19.2,22.4,25.4,28.3,33.7,38.0,41.2,43.1,44.2]
bolsig01=[]
bolsig1=[]
bolsig5=[]
for j in range(len(ArRat)):
    REr=[]
    haddad=[]
    if j==0:
        REr=REr01
        haddad=haddad01
        vDriftFile = open('vdriftBia2Ar999N001.out','w')
    elif j==1:
        REr=REr1
        haddad=haddad1
        vDriftFile = open('vdriftBia2Ar990N010.out','w')
    elif j==2:
        REr=REr5
        haddad=haddad5
        vDriftFile = open('vdriftBia2Ar950N050.out','w')
     
    for i in range(len(REr)):
        InputMaker(GasTemperature,REr[i],ArRat[j],NRat[j])
        Runner()
        electronRawMobility=OutputExtractor()
        driftVelocity=electronRawMobility*REr[i]*1e-21
        vDriftFile.write(str(REr[i])+'\t'+str(driftVelocity)+'\t'+str(haddad[i]*1e3)+'\n')
        if j==0: bolsig01.append(driftVelocity/1e3)
        elif j==1: bolsig1.append(driftVelocity/1e3)
        elif j==2: bolsig5.append(driftVelocity/1e3)
    vDriftFile.close()    

import matplotlib.pyplot as plt
plt.plot(REr01,haddad01,'r.')
plt.hold(True)
plt.plot(REr1,haddad1,'g.')
plt.plot(REr5,haddad5,'b.')
plt.plot(REr01,bolsig01,'r-')
plt.plot(REr1,bolsig1,'g-')
plt.plot(REr5,bolsig5,'b-')
plt.xlabel('E/N (Td)')
plt.ylabel('Drift velocity (m/s)')
plt.savefig('gasdriftveloArRatW.png', format='png', dpi=1000)
