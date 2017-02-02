# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 2017

@author: Andoni Torres
"""

#Programme qui lit les données de graphiques lsd et délimite les bornes d'intégration selon le fit gaussien de Stokes I et calcule l'intrégrale entre ces bornes des paramètres de Stokes afin de calculer le champ longitudinal Bl et Nl

import numpy as np
import os 
import glob
from scipy.integrate import simps
from scipy.optimize import curve_fit

os.chdir("/data/coolsnap/data_coolsnap_archive")
#os.chdir("/data/coolsnap/spectra_espadons_coolsnap/")

def gaus(x,a,x0,sigma,Ic):
    return Ic-a*np.exp(-(x-x0)**2/(2*sigma**2))
    
def integ(a,b):
    return float(simps(a,b))

liste = []
atrocites = [] # Fichiers dont le graphique Stokes I ne représentent pas une gausienne
compteur = 0

for file in glob.glob("[0-9]*pn.lsd"):
   liste.append(file)

wave_lande = np.genfromtxt('/h/torres/Downloads/mean_wave_lande(true)',skip_header=0, dtype=[('0','S15'),('1','f8'),('2','f8')])

noms = wave_lande['0'] 
wave = wave_lande['1']
lande = wave_lande['2']

"""# Lire les noms des fichiers a partir d'une liste
liste = np.genfromtxt('/h/torres/toto',dtype=str)"""

print ('Fichier', 'Bl', 'erreur', 'x0',	'sigma', 'Ic',	'Hauteur du pic')

for filename in liste:

    for i in range(0,(len(noms))):
        if filename == noms[i]:
            geff = lande[i]
            lambda0 = wave[i]
        else:
            geff = 1.232
            lambda0 = 668.1893
	cc = 299792
	coeff = (-2.14*10**11)/(lambda0*cc*geff)

#    os.chdir("/data/coolsnap/data_coolsnap_archive")
    A = np.genfromtxt(filename,skip_header=2, dtype=[('0','f8'),('1','f8'),('2','f8'),('3','f8'),('4','f8'),('5','f8'),('6','f8')])
    
# Defnir les donnes
# Vitesse radiale
    VR = A['0']
# Stokes I
    SI = A['1']
# Stokes V
    SV = A['3']
# Stokes N
    SN = A['5']
    
# Erreurs sur les données
    eSI = A['2']
    eSV = A['4']
    eSN = A['6']
       
# Lissage gaussien (permet de définir les bornes d'intégration)
      
    Xmoy = np.mean(VR[np.where(SI == min(SI))]) # Abcisse du SI minimum
    
    params0=[0.2, Xmoy, 5.0, 1.0]  # Estimation initiale des paramètres
    
    popt, pcov = curve_fit(gaus, VR, SI, params0) # popt[]:paramètres optimisés; popt[0]: ordonnée moyenne ,popt[1]: abscisse moyenne, popt[2]:sigma, popt[3]: Icontinu
    
    """#Graphique
    
    x = np.arange(-243,243, 0.2)
    y = gaus(x, popt[0], popt[1], popt[2], popt[3])

    import matplotlib.pyplot as plt
    plt.plot(VR,SI,'k.',label='donnees')
    plt.plot(x,y,'r-',label='lissage')
    plt.show()"""
            
# Borne inferieure et supérieure d'intégration (+/- 3 sigma)
    binf = popt[1]-3*popt[2]
    bsup = popt[1]+3*popt[2]
    
    if binf > -243 and bsup < 243:
# Créer un array des index du domaine d'intégration 
        bornesVR = ([np.where((VR >= binf) & (VR <= bsup))])
        domaine = np.array(bornesVR)

# Données à intégrer
        bVR = VR[domaine]
        bSI = popt[3]-SI[domaine] #Ic-I(vr)
        bSV = bVR*SV[domaine] #vr*V(vr)
        bSN = SN[domaine]
            
# Intégration
        iSI = integ(bSI,bVR)
        iSV = integ(bSV,bVR)
        iSN = integ(bSN,bVR)
        
# Champ longitudinal
        Bl = float(coeff*(iSV/iSI))

# Calcul de l'erreur
        ebSI = eSI[domaine]
        ebSV = eSV[domaine] 
        ebSN = eSN[domaine]

        eiSI = integ(ebSI,bVR)
        eiSV = integ(ebSV,bVR)
        eiSN = integ(ebSN,bVR)
        
        erSV = ((coeff/iSI)*(eiSV**2))**2 
        erSI = ((coeff*iSV/(iSI**2))*(eiSI**2))**2
        erreur = float(np.sqrt(erSV+erSI))

    	Xmin = np.mean(VR[np.where(SV == min(SV[domaine]))]) # Abcisse du SV minimum
		print(Xmin)

        compteur = compteur +1

#        print(filename, "%.f" % Bl, "%.f" %  erreur, "%.f" %  popt[1],	"%.4f" %  popt[2],	"%.4f" %  popt[3],	"%.4f" %  (popt[3]-popt[0]))
#        print filename, '\t', "%.f" % Bl,'\t', "%.f" %  erreur,'\t', "%.f" %  popt[1],'\t',	"%.4f" %  popt[2],'\t',	"%.4f" %  popt[3],'\t',	"%.4f" %  (popt[3]-popt[0])
    else:
		atrocites.append(filename)

#print(len(liste), compteur, atrocites)
