# -*- coding: utf-8 -*-
"""Created on Tue Jan 17 2017 @author: Andoni Torres Programme qui lit les données de graphiques lsd et délimite les
bornes d'intégration selon le fit gaussien de Stokes I et calcule l'intrégrale entre ces bornes des paramètres de
Stokes afin de calculer le champ longitudinal Bl et Nl """

import numpy as np
import os
import glob
from scipy import interpolate
from scipy.integrate import simps
from scipy.optimize import curve_fit


def gaus(x, a, x0, sigma, ic):
    return ic - a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def integ(a, b):
    return float(simps(a, b))

def func(x, a, x0, sigma, ic, a2, x02, sigma2, ic2):
    return gaus + gaus
params01 = [0.1, -0.1, 0.01, 0]
params02 = [-0.1, 0.1, 0.01, 0]

liste = []
atrocites = []  # Liste des fichiers dont le graphique Stokes I ne représentent pas une gausienne
cc = 299792

for file in glob.glob("/data/coolsnap/data_coolsnap_archive/[0-9]*pn.lsd"):
#for file in glob.glob("/data/coolsnap/spectra_espadons_coolsnap/[0-9]*pn.lsd"):
    liste.append(file)

wave_lande = np.genfromtxt('/h/torres/Downloads/mean_wave_lande(true)', dtype=[('0', 'S15'), ('1', 'f8'), ('2', 'f8')])

noms = wave_lande['0']
wave = wave_lande['1']
lande = wave_lande['2']

"""# Lire les noms des fichiers a partir d'une liste
liste = np.genfromtxt('/h/torres/toto',dtype=str)"""

print 'Champ longitudinal: Fichier', 'Bl', 'erreur', 'Paramètres de Stokes V: min(x,y)', 'min(x,y)'

for filename in liste:

    for i in range(0, (len(noms))):
        if filename == noms[i]:
            geff = lande[i]
            lambda0 = wave[i]
        else:
            geff = 1.232
            lambda0 = 668.1893
        coeff = (-2.14 * 10 ** 11) / (lambda0 * cc * geff)

    #    os.chdir("/data/coolsnap/data_coolsnap_archive")
    A = np.genfromtxt(filename, skip_header=2,
                      dtype=[('0', 'f8'), ('1', 'f8'), ('2', 'f8'), ('3', 'f8'), ('4', 'f8'), ('5', 'f8'), ('6', 'f8')])

    # Defnir les donnes

    VR = A['0']  # Vitesse radiale
    SI = A['1']  # Stokes I
    SV = A['3']  # Stokes V
    SN = A['5']  # Stokes N

    # Erreurs sur les données
    eSI = A['2']
    eSV = A['4']
    eSN = A['6']

    # Lissage gaussien (permet de définir les bornes d'intégration)

    Xmoy = np.mean(VR[np.where(SI == min(SI))])  # Abcisse du SI minimum

    params0 = [0.2, Xmoy, 5.0, 1.0]  # Estimation initiale des paramètres de la gaussienne

    popt, pcov = curve_fit(gaus, VR, SI,
                           params0)  # popt[]:paramètres optimisés; popt[0]: ordonnée moyenne ,popt[1]: abscisse
    # moyenne, popt[2]: écart-type, popt[3]: continuum

    """#Graphique
    
    x = np.arange(-243,243, 0.2)
    y = gaus(x, popt[0], popt[1], popt[2], popt[3])

    import matplotlib.pyplot as plt
    plt.plot(VR,SI,'k.',label='donnees')
    plt.plot(x,y,'r-',label='lissage')
    plt.show()"""

    # Borne inferieure et supérieure d'intégration (+/- 3 sigma)
    binf = popt[1] - 3 * popt[2]
    bsup = popt[1] + 3 * popt[2]

    # Création d'un array des index du domaine d'intégration
    if binf > -243 and bsup < 243:
        bornesVR = ([np.where((VR >= binf) & (VR <= bsup))])
        domaine = np.array(bornesVR)

        # Données à intégrer
        bVR = VR[domaine]
        bSI = popt[3] - SI[domaine]  # continuum - SI(VR)
        bSV = bVR * SV[domaine]  # VR*V(VR)
        bSN = SN[domaine]

        # Intégration (par la loi de Simpson)
        iSI = integ(bSI, bVR)
        iSV = integ(bSV, bVR)
        iSN = integ(bSN, bVR)

        # Champ longitudinal
        Bl = float(coeff * (iSV / iSI))


        if len(SV[domaine][0][0]) > 2:
            k=3
        else:
            k=1

        s = interpolate.UnivariateSpline(bVR[0][0], SV[domaine][0][0], k=k, s=0)

        xs = np.arange(bVR[0][0][0], bVR[0][0][-1], 0.01)
        y = s(xs)

        abc = max(y)
        compton = float(xs[np.where(y == max(y))])
        efg = min(y)
        bompton = float(xs[np.where(y == min(y))])

        import matplotlib.pyplot as plt

        """plt.figure()
        plt.plot(bVR[0][0], SV[domaine][0][0], 'k.', label='donnees')
        plt.plot(xs, y, 'r-', label='lissage')
        plt.show()"""

        # Calcul de l'erreur
        ebSI = eSI[domaine]
        ebSV = eSV[domaine]

        eiSI = integ(ebSI, bVR)
        eiSV = integ(ebSV, bVR)

        erSV = ((coeff / iSI) * (eiSV ** 2)) ** 2
        erSI = ((coeff * iSV / (iSI ** 2)) * (eiSI ** 2)) ** 2
        erreur = float(np.sqrt(erSV + erSI)) + abs(iSN)
    	print filename, '\t', "%.f" % Bl,'\t', "%.f" %  erreur,'\t',"%.5f" % abc,'\t',"%.2f" % compton,'\t',"%.5f" % efg,'\t',"%.2f" % bompton

    #        print(filename, "%.f" % Bl, "%.f" %  erreur, "%.f" %  popt[1],	"%.4f" %  popt[2],	"%.4f" %  popt[3],	"%.4f" %  (popt[3]-popt[0]))
    #        print filename, '\t', "%.f" % Bl,'\t', "%.f" %  erreur,'\t', "%.f" %  popt[1],'\t',	"%.4f" %  popt[2],'\t',	"%.4f" %  popt[3],'\t',	"%.4f" %  (popt[3]-popt[0])
    else:
        atrocites.append(filename)

# print(len(liste), atrocites)
