# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 2017

@author: Andoni Torres
"""

#Programme qui lit les données de certains graphiques lsd et délimite les bornes d'intégration selon le fit gaussien de Stokes I et calcule l'intrégrale entre ces bornes des graphiques des paramèetres de Stokes afin de calculer le champ longitudinal Bl et Nl et les retourner en une liste
#Peut-etre rajouter commandes pour comparer avec programme de Claire et sortir l'ecart

import numpy as np
import os 
import glob
from scipy.integrate import simps
from scipy.optimize import curve_fit

lambda0 = 750
geff = 1.25
c = 299792

liste = []
os.chdir("C:/Users/Andoni Torres/Desktop/Stage_Hiver_2017/Champ_longitudinal")
for file in glob.glob("[0-9]*.lsd"):
    liste.append(file)
    
#liste = np.genfromtxt('noms.txt',dtype=str)

for filename in liste:
    os.chdir("C:/Users/Andoni Torres/Desktop/Stage_Hiver_2017/Champ_longitudinal")
    A = np.genfromtxt(filename,skip_header=2, usecols=(0,1,3,5),dtype=[('0','f8'),('1','f8'),('3','f8'),('5','f8')])
    
#Defnir les donnes
#vitesse radiale
    vr = A['0']
#Stokes I
    sI = A['1']
#Stokes V
    sV = A['3']
#Stokes N
    sN = A['5']
    
    min1= min(sI)
    min2=np.where(sI == min1)
    min3=vr[min2]
       
#Lissage gaussien    
    
    def gaus(x,a,x0,sigma,Ic):
        return Ic-a*np.exp(-(x-x0)**2/(2*sigma**2))
        
    params0=[0.1,0,min3,1.0]  # Estimation initiale des paramètres
    

    popt, pcov = curve_fit(gaus, vr, sI, params0) # popt[]:paramètres optimisés; popt[0]: ordonnée moyenne ,popt[1]: abscisse moyenne, popt[2]:sigma, popt[3]: Icontinu
        
    x = np.arange(-243,243, 0.2)
    y = gaus(x, popt[0], popt[1], popt[2], popt[3])
    
# Graphique (non nécessaire)
    import matplotlib.pyplot as plt
    plt.plot(vr,sI,'k.',label='donnees')
    plt.plot(x,y,'r-',label='lissage')
    plt.show()
            
# Borne inferieure et supérieure d'intégration (+/- 3 sigma)
    binf=popt[1]-3*abs(popt[2])
    bsup=popt[1]+3*abs(popt[2])
    
    if binf>-243 and bsup<243:
# Creer un array des index du domaine d'intégration 
        bornesVR=([np.where((vr>=binf) & (vr<=bsup))]) 
        
# Données à intégrer
        bVR=vr[(np.array(bornesVR))]
        bsI=popt[3]-sI[(np.array(bornesVR))]
        bsV=vr[(np.array(bornesVR))]*sV[(np.array(bornesVR))]
        bsN=sN[(np.array(bornesVR))]
            
# Calcul des intégrales  
        Ii=float(simps(bsI,bVR))
        Iv=float(simps(bsV,bVR))
        In=float(simps(bsN,bVR))
        
# Champ longitudinal
        Bl= float((-2.14*10**11*(Iv/Ii))/(lambda0*c*geff))
        
        print(filename, Bl)
    else:
        pass
