# -*- coding: utf-8 -*-
# matplotline inline
"""
Created on Tue Jan 17 2017 @author: Andoni Torres

Programme qui lit les données de graphiques lsd et délimite les bornes d'intégration selon le fit gaussien de Stokes I
et calcule l'intrégrale entre ces bornes des paramètres de Stokes afin de calculer le champ longitudinal Bl
"""

import numpy as np
from scipy import interpolate
from scipy.integrate import simps
from scipy.optimize import curve_fit
import Dictionnaire
import matplotlib.pyplot as plt
import lmfit

#plt.ion()

Etoiles = Dictionnaire.etoiles() # Dictionnaire des étoiles et de leurs paramètres
Stokes = Dictionnaire.stokes() # Dicitonnaire des odomètres et des paramètres des profils Stokes

# Note: Changer pour une gausienne + lorentzienne
def gaus(x, a, x0, sigma, ic):
    return ic - a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

def integ(a, b): # Intégration point-par-point par la loi de Simpson
    return float(simps(a, b))


liste = []
atrocites = []  # Liste des fichiers dont le graphique Stokes I ne représentent pas une gausienne
cc = 299792  # Vitesse de la lumière km/s

print 'Champ longitudinal: Fichier', 'Bl', 'erreur', 'Paramètres de Stokes V: min(x,y)', 'min(x,y)'

for star, param in Etoiles.items():
    for item in param[0]:
        if item in Stokes.keys():
            VR = np.array(Stokes[item][0]) # Vitesse radiale
            SI = np.array(Stokes[item][1]) # Stokes I
            SV = np.array(Stokes[item][2]) # Stokes V
            SN = np.array(Stokes[item][3]) # Stokes N

            eSI = np.array(Stokes[item][4]) # Erreurs sur les données
            eSV = np.array(Stokes[item][5])
            eSN = np.array(Stokes[item][6])

            lambda0 = Stokes[item][7] # Facteur de landé
            geff = Stokes[item][8] # g effectif

            coeff = (-2.14 * 10 ** 11) / (lambda0 * cc * geff)

            # Lissage gaussien (permet de définir les bornes d'intégration)

            if SI.tolist(): # Si SI n'est pas vide
                #Xmoy = np.mean(VR[SI.index(min(SI))])  # Abcisse du SI minimum

                Xmoy = np.mean(VR[np.where(SI == min(SI))])  # Abcisse du SI minimum

                params0 = [0.2, Xmoy, 5.0, 1.0]  # Estimation initiale des paramètres de la gaussienne

                # popt[]:paramètres optimisés; popt[0]: ordonnée moyenne ,popt[1]: abscisse moyenne, popt[2]: écart-type, popt[3]: continuum
                popt, pcov = curve_fit(gaus, VR, SI, params0)

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
                        k = 3
                    else:
                        k = 1

                    s = interpolate.UnivariateSpline(bVR[0][0], SV[domaine][0][0], k=k, s=0)
                    xs = np.arange(bVR[0][0][0], bVR[0][0][-1], 0.01)
                    ys = s(xs)

                    abc = max(ys) # max de l'interpolation
                    toto1 = float(xs[np.where(ys == max(ys))]) # abscisse max
                    efg = min(ys) # min de l'interpolation
                    toto2 = float(xs[np.where(ys == min(ys))]) # abscisse min

                    # Calcul de l'erreur
                    ebSI = eSI[domaine]
                    ebSV = eSV[domaine]

                    eiSI = integ(ebSI, bVR)
                    eiSV = integ(ebSV, bVR)

                    erSV = ((coeff / iSI) * (eiSV ** 2)) ** 2
                    erSI = ((coeff * iSV / (iSI ** 2)) * (eiSI ** 2)) ** 2
                    erreur = float(np.sqrt(erSV + erSI)) + abs(iSN)




                    if len(param[0]) % 2 == 0:
                        col, rang = len(param[0])/2, 2
                    elif len(param[0]) == 1:
                        col, rang = 1, 1
                    else:
                        col, rang = len(param[0])/2+1, 2

                    ax = plt.subplot(col, rang, param[0].index(item) + 1)

                    x = np.arange(-243, 243, 0.2)
                    y = gaus(x, popt[0], popt[1], popt[2], popt[3])

                    ax.plot(VR,SI,'k.') # Stokes I sur tout le domaine
                    #ax.plot(bVR, SI[domaine], 'r-') # Gausienne

                    ax.fill_between(VR, SI, popt[3], where=binf < VR, facecolor='green', interpolate=True)
                    #ax.plot(bVR[0][0], SV[domaine][0][0], 'k.') # Stokes V sur le domaine d'intégration
                    #ax.plot(xs, ys, 'r-') # Interpolation de Stokes V

                    ax.set_title(str(item)+', '+ str(param[5][param[0].index(item)]))

                    #ax.plot([toto1, toto1], [0, abc], color='k', linestyle='--')
                    #ax.text(0, abc, str("aaaa"), fontsize=16)

                    #ax.plot([toto2, toto2], [0, efg], color='k', linestyle='--')

                    #ax.axhline(y=0, color='k', linestyle='-')
                    #ax.axvline(x=0, color='k', linestyle='-')



                    print item, '\t', "%.2f" % Bl, '\t', "%.2f" % erreur, '\t', "%.5f" % abc, '\t', "%.2f" % toto1, '\t', "%.5f" % efg, '\t', "%.2f" % toto2
                else:
                    atrocites.append(item)
            else:
                atrocites.append(item)
    #plt.tight_layout()
    plt.suptitle('Stokes N and interpolation curve for object ' + str(star), fontsize=16)
    plt.show()

def main():
    pass

if __name__ == "__main__":
    main()
