# -*- coding: utf-8 -*-

"""
Created on Tue Feb 15 2017 @author: Andoni Torres

"""

import numpy as np
import Dictionnaire

A = np.genfromtxt('/data/coolsnap/tables/Zeembroad.csv', delimiter=',', usecols=(1, 4),
                  dtype=[('0', 'S15'), ('1', 'f8')])

odo = A['0']
zeem = A['1']
odo_upd = []

for nom in odo:
    odo_u = nom.replace("p","").replace("i", "")
    odo_u = odo_u.replace(" ", "")
    odo_upd.append(odo_u)

Etoiles = Dictionnaire.etoiles()

noms = [] # liste des odomètres dans le dicitonnaire Etoiles
zams = [] # liste des facteurs de zeeman
Bl = [] # Liste des Bl associés aux odomètres
err = [] # Erreurs sur Bl

compteur = 0

for i, v in Etoiles.items():
    for w in range (0, len(v[0]), 1):
        if v[0][w] in odo_upd:
            if v[1][w] != 'nulle':
                ind = odo_upd.index(v[0][w])
                noms.append(odo_upd[ind])
                zams.append(zeem[ind])
                Bl.append(abs(v[3][w]))
                err.append(v[4][w])
                compteur = compteur + 1

import matplotlib.pyplot as plt

xerr = 0.5
yerr = err
fig, ax = plt.subplots()

plt.xlim([0, 7.5])
plt.ylim([0, 50])
plt.xlabel('Zeeman broadening')
plt.ylabel('Champ longitudinal absolu (G)')
ax.grid(True, which='both')
ax.axhline(y=0, color='k', linestyle='-')
ax.axvline(x=0, color='k', linestyle='-')
ax.errorbar(zams, Bl, xerr=xerr, yerr=yerr, color='k', fmt='.')
plt.show()