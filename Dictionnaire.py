# -*- coding: utf-8 -*-

"""
Ce programme crée deux dictionnaires.
'Etoiles' est un dictionnaire avec le nom des objets stellaires en tant que Key et qui
a comme Value un array avec des listes de paramètres
'Stokes' a pour Key les odomètres et comme Value un array avec des listes des données des
profils de Stokes et leurs incertitudes.
"""

import numpy as np
import os
import glob


def etoiles():
    liste_coolsnap = np.genfromtxt('/h/www/www.cfht.hawaii.edu/coolsnap/clichesfroids_log.dat',
                                   delimiter=',', usecols=(3, 2, 16, 24),
                                   dtype=[('0', 'S25'), ('1', 'S25'), ('2', 'S25'), ('3', 'f8')])

    liste_archive = (np.genfromtxt('/h/www/www.cfht.hawaii.edu/coolsnap/clichesfroids_log_archive_p.dat',
                                   delimiter=',', usecols=(3, 2, 17, 31),
                                   dtype=[('0', 'S25'), ('1', 'S25'), ('2', 'S25'), ('3', 'f8')]))

    objet = np.concatenate((liste_archive['0'], liste_coolsnap['0']), axis=0)
    odometre = np.concatenate((liste_archive['1'], liste_coolsnap['1']), axis=0)
    detection = np.concatenate((liste_archive['2'], liste_coolsnap['2']), axis=0)
    vsini = np.concatenate((liste_archive['3'], liste_coolsnap['3']), axis=0)

    # Crée une nouvelle liste des noms des étoiles standarisé sans espaces et en majuscules et odomètres sans espaces

    liste_update = []
    global odo_upd
    odo_upd = []

    for i in objet:
        star = i.replace(" ", "")
        star = star.upper()
        liste_update.append(star)

    for i in odometre:
        odu = i.replace(" ", "")
        odo_upd.append(odu)

    # Crée un dictionnaire qui associe chaque étoile à un array de paramètres
    # Etoiles[nom de l'étoile] = [[Odometre], [Detection], [vsini], [Bl]]

    Etoiles = {}

    for i in range(0, (len(liste_update))):
        Etoiles[liste_update[i]] = [[],[],[],[],[]]

    for star in Etoiles.keys():
        for i in range(0, (len(liste_update))):
            if star == liste_update[i]:
                Etoiles[star][0].append(odometre[i].replace(" ", ""))
                Etoiles[star][1].append(detection[i].replace(" ", ""))
                Etoiles[star][2].append(vsini[i])

    # PARTIE BL

    liste_bl = np.genfromtxt('/h/torres/Downloads/Bl.txt', usecols=(0, 1, 2), skip_header=1,
                             dtype=[('0', 'S25'), ('1', 'f8'), ('2', 'f8')])

    filenames = liste_bl['0']
    Bl = liste_bl['1']
    blerr = liste_bl['2']


    for i , v in Etoiles.items():
        for w in v[0]:
            if w in filenames:
                ind = np.where(filenames == w)
                Etoiles[i][3].append(Bl[ind[0][0]])
                Etoiles[i][4].append(blerr[ind[0][0]])
            else:
                Etoiles[i][3].append(np.nan)
                Etoiles[i][4].append(np.nan)

    return Etoiles


# Crée un dictionnaire qui associe chaque odomètre à un array de cordonnées de profils Stokes
# Stokes[odomètre] = [[VR], [Stokes I], [Stokes V], [Stokes N],
#                     [Incertitude Stokes I], [Incertitude Stokes V], [IncertitudeStokes N],
#                     facteur de landé, g effectif]


def stokes():
    #668.1893
    Stokes = dict.fromkeys(odo_upd, [[],[],[],[],[],[],[], np.nan, np.nan])

    paths = ["/data/coolsnap/data_coolsnap_archive/", "/data/coolsnap/spectra_espadons_coolsnap/"]
    liste_odo = []

    for path in paths:
        os.chdir(path)
        for file in glob.glob("[0-9]*pn.lsd"):
            if file.replace("pn.lsd", "") in Stokes.keys():

                A = np.genfromtxt(file, skip_header=2,
                                  dtype=[('0', 'f8'), ('1', 'f8'), ('2', 'f8'), ('3', 'f8'), ('4', 'f8'), ('5', 'f8'),
                                         ('6', 'f8')])

                Stokes[file.replace("pn.lsd", "")] = [A['0'].tolist(),A['1'].tolist(),A['2'].tolist(),A['3'].tolist(),
                                                      A['4'].tolist(),A['5'].tolist(),A['6'].tolist(),668.1893, 1.232]

    # Landé et geff

    wave_lande = np.genfromtxt('/h/torres/Downloads/mean_wave_lande(true)',
                               dtype=[('0', 'S15'), ('1', 'f8'), ('2', 'f8')])

    noms = wave_lande['0']
    wave = wave_lande['1']
    lande = wave_lande['2']

    for i in range(0, len(noms), 1):
        Stokes[noms[i].replace("pn.lsd", "")][7] = wave[i]
        Stokes[noms[i].replace("pn.lsd", "")][8] = lande[i]

    return Stokes

def main():
    pass

if __name__ == "__main__":
    main()
