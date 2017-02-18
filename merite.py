# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 2017

@author: Andoni Torres
"""

# Ébauche de programme pour la fonction de mérite

# Accorder une pondération à chaque paramètres de l'étoile

import numpy as np
import math
import Dictionnaire

Etoiles = Dictionnaire.etoiles() # Importer les donnés

# Donner des points pour les paramètres
# 10 points au total, vsini: 3 points, Champ longitudinal: 3 points, Detections: 4 points

ponderation = [1, 8, 1]

# Génération du coefficient du paramètre vsini

liste1 = []  # Liste des moyennes de vsini pour histogramme
liste2 = []  # Liste de tous les coefficients pour histogramme

coefficients = dict.fromkeys(Etoiles, [])


for i, v in Etoiles.items():
    liste_non_nan = [valeur for valeur in v[2] if not math.isnan(valeur)]  # Liste des vsini sans les valeurs 'nan'
    if liste_non_nan != []:
        moyenne = np.mean(liste_non_nan)
        liste1.append(moyenne)
        if 0 < moyenne < 2:
            coeff = 1.0
        elif 2 < moyenne < 2.5:
            coeff = 1.0 - 0.8
        elif 2.5 < moyenne < 3:
            coeff = 1.0 - 0.75 / 4
        elif 3 < moyenne < 4:
            coeff = 1.0 - 0.75 / 3
        elif 5 < moyenne < 6:
            coeff = 1.0 - 0.75 / 2
        liste2.append(coeff)
        coefficients[i] = [coeff]
    else:
        coeff = 0.15
        liste2.append(coeff)
        coefficients[i] = [coeff]

"""
# Histogramme des vsini moyens (valeur absolue) et des coefficents

import matplotlib.pyplot as plt

binf = -max(liste1)
bsup = max(liste1)

bins = np.arange(0, bsup, 1)

plt.xlim([binf, bsup])
plt.hist(liste1, bins=bins, alpha=0.5)
plt.show()

plt.xlim([0, 1])
plt.hist(liste2, alpha=0.5)
plt.show()
"""

# Génération du coefficient du paramètre Bl

intensite = []
liste3 = []

for i, v in Etoiles.items():
    for item in v[3]:
        if len(v[3]) > 1:
            valeur1 = abs(max(v[3]) - min(v[3]))
            intensite.append(valeur1)
        else:
            valeur1 = abs(v[3][0])
            intensite.append(valeur1)

        if 0 < valeur1 < 15:
            coeff2 = 1.0
        elif 15 < valeur1 < 70:
            coeff2 = 0.7
        elif 70 < valeur1 < 200:
            coeff2 = 0.4
        elif 200 < valeur1 < 500:
            coeff2 = 0.2
        elif 500 < valeur1:
            coeff2 = 0
    liste3.append(coeff2)
    coefficients[i].append(coeff2)

"""
# Histogramme des Champ longitudinaux (Bl) et du coefficient

import matplotlib.pyplot as plt

binf = 0
bsup = 2000

#bins = np.arange(0, bsup, 25)
# popt, pcov = curve_fit(gaus, abcs, plt.hist(liste4, bins=bins, alpha=0.5), params0)

#plt.xlim([binf, bsup])
#plt.hist(intensite, bins=bins, alpha=0.5)
#plt.show()

plt.xlim([0, 1])
plt.hist(liste3, alpha=0.5)
plt.show()
"""

# Génération du coefficient du paramètre Detection
# Détection Definite = 3 points, Marginale = 2 points, Nulle = -2 points

liste4 = []

for i, v in Etoiles.items():
    points = 0
    for detect in v[1]:
        if detect == 'nulle':
            points = points + 3.0
        elif detect == 'marginale':
            points = points - 2.0
        elif detect == 'definite':
            points = points - 3.0
            if points < 0:
                points = 0
    coeff3 = points / (3 * len(v[1]))
    liste4.append(coeff3)
    coefficients[i].append(coeff3)

"""
# Histogramme du coefficient pour les détections

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.xlim([0, 1])
plt.hist(liste4, alpha=0.5)
plt.show()
"""

# Accorder une cote à chaque étoile dépendament de ses paramètres et des pondérations

etoiles_cotes = dict.fromkeys(Etoiles, )
liste_cotes = []

for i,v in coefficients.items():
    cote = sum(np.multiply(v, ponderation))
    liste_cotes.append(cote)
    etoiles_cotes[i] = cote
    print i, v


for i in sorted(etoiles_cotes.items(), key=lambda x:x[1], reverse=True):
    print i[0],':', "%.2f" % i[1]


# Histogramme des cotes

import matplotlib.pyplot as plt

plt.xlim([0, 10])
plt.hist(liste_cotes, alpha=0.5)
plt.show()

def main():
    pass

if __name__ == '__main__':
  main()