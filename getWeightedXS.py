# -*- coding: utf-8 -*-
"""
Quick script too lookup fission cross section data, weight it with an analytic 
spectrum, and produce a similar figure as in NIMA Vol 593 pp510-518

Zs. Elter, 2017. march 20.
"""


import urllib.request
import matplotlib.pyplot as plt
import numpy as np

def getXS(content):
    """getting data from the html content. I tried to avoid Beautiful soup.
    the content has a shape as follows like:    
    Energy(eV) XS(b)<br>
    1.00000E-05 3.07139<br>"""
    energy=[]
    xs=[]
    flag=False
    i=0
    for line in content:
        x=line.strip().split()
        if  x[0]==b'Energy(eV)':
            flag=True
            continue
        if x[0]==b'</span>':
            flag=False
        if flag:
            energy.append(float(x[0]))
            xs.append(float(x[1][:-4]))
    return energy, xs

def getSpectrum(energy,xs):
    """Approximating the spectrum as in NIMA Vol 593 pp510-518 (Philippe's) paper. 
    Constants describing the flux are taken from the same paper"""
    fth=[]
    k=8.6173324e-5 #eV/K
    T=325 #K
    for E in energy:
        fth.append(E*np.exp(-E/(k*T)))   #who needs pythonic way?
    energy=np.array(energy)  #i want np arrays, though I don't like np.append.
    fth=np.array(fth)
    fth=fth/np.trapz(fth,energy)
    
    ##epithermal
    fepi=[]
    E0=0
    E1=0.2
    E2=0.25e6
    E3=1e6
    for E in energy:
        if E0<=E and E<E1:
            fepi.append((E**2-E0**2)/(E*(E1**2-E0**2)))
        elif E1<=E and E<E2:
            fepi.append(1/E)
        elif E2<=E and E<E3:
            fepi.append((E**2-E3**2)/(E*(E2**2-E3**2)))
        else:
            fepi.append(0.0)
    fepi=np.array(fepi)
    fepi=fepi/np.trapz(fepi,energy)
    
    #fast
    ffast=[]
    a=8.09e5 #u235 9.65e5
    b=9.32e-7 #u235 2.29e-6
    for E in energy:
        ffast.append(np.exp(-E/a)*np.sinh(np.sqrt(b*E)))
    ffast=np.array(ffast)
    ffast=ffast/np.trapz(ffast,energy)
    
    #create real spectra
    PhiTh=3.2e14
    PhiEpi=1.98e14
    PhiFast=0.98e14
    
    phi=PhiTh*fth+PhiEpi*fepi+PhiFast*ffast
    
    #sigmabar
    xs=np.array(xs)
    sigmabarF=np.trapz(xs*(PhiFast*ffast),energy)/np.trapz((PhiFast*ffast),energy)
    sigmabarTh=np.trapz(xs*(PhiTh*fth),energy)/np.trapz((PhiTh*fth),energy)

    sigmabarFTh=sigmabarF/sigmabarTh
    return sigmabarF,sigmabarFTh

#####main script

#           Name, Z,    A,   MAT
isotopes=[['Th', '90', '232','9040'],
          ['U', '92', '238','9237'],
          ['U', '92', '236','9231'],
          ['Pa', '91', '231','9131'],
          ['Pu', '94', '242','9446'],
          ['Np', '93', '237','9346'],
          ['Am', '95', '243','9549'],
          ['Cm', '96', '246','9643'],
          ['Pu', '94', '240','9440'],
          ['Cm', '96', '248','9649'],
          ['U', '92', '234','9225'],
          ['Am', '95', '241','9543'],
          ['Cm', '96', '244','9637'],
          ['U', '92', '235','9228'],
          ['Cf', '98', '252','9861'],
          ['Pu', '94', '238','9434'],
          ['Pu', '94', '241','9443'],
          ['Pu', '94', '239','9437'],
          ['U', '92', '233','9222'],
          ['Cm', '96', '247','9646'],
          ['Pu', '94', '236','9428'],
          ['U', '92', '232','9219'],
          ['Cm', '96', '243','9634'],
          ['Cf', '98', '251','9858'],
          ['Cm', '96', '245','9640']]

website='http://atom.kaeri.re.kr/nuchart/getData.jsp?target=jeff3.2,%s,%s,%s,3,18'

sigmaf=[]
sigmafth=[]
import time
for i in range(len(isotopes)):    
    url=website%(isotopes[i][1],isotopes[i][2],isotopes[i][3])
    if i==0:  #Th232 not included in jeff 3.2
        url='http://atom.kaeri.re.kr/nuchart/getData.jsp?target=jendl4.0,90,232,9040,3,18'
    print('retrieving '+url+' ...')
    with urllib.request.urlopen(url) as response:
        content = response.readlines()
    energy,xs=getXS(content)
    #spectrum is computed each time because the energy binning is not necessarily the same
    sigmafi,sigmafthi=getSpectrum(energy,xs)
    sigmaf.append(sigmafi)
    sigmafth.append(sigmafthi)
    time.sleep(10) #not to overwhelm the server:)


###plot
fig, ax = plt.subplots()
ax.semilogy(sigmaf,sigmafth,'x')
ax.hold(True)
ax.semilogy(sigmaf[1],sigmafth[1],'r.',markersize=15)
ax.semilogy(sigmaf[4],sigmafth[4],'r.',markersize=15)
ax.semilogy(sigmaf[13],sigmafth[13],'r.',markersize=15)
ax.semilogy(np.linspace(0,3),np.ones((len(np.linspace(0,3)),1)),'k--')
    
for j in range(len(isotopes)):
    ax.annotate('  '+isotopes[j][0]+isotopes[j][2], (sigmaf[j],sigmafth[j]))
ax.set_xlabel(r'$\bar{\sigma}_f$ (barns)')
ax.set_ylabel(r'$\bar{\sigma}_f/\bar{\sigma}_{th}$')
ax.set_xlim([0,2.75])
ax.set_ylim([1e-4,3e4])
plt.savefig('FCxsscatter.png', format='png', dpi=1000)
