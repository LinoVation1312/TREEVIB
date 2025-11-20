# -*- coding: utf-8 -*-
"""
Script : tracer plusieurs variables météo depuis un fichier Excel (Corrigé)
"""

import pandas as pd
import matplotlib.pyplot as plt
from tkinter import Tk, filedialog
import numpy as np # Importation de numpy pour une meilleure gestion des types et des erreurs

# --- 1. Sélection des fichiers Excel ---
root = Tk()
root.withdraw() # cache la fenêtre principale Tkinter

# On demande le premier fichier
fichier_excel_meteostat = filedialog.askopenfilename(
    title="Choisir le fichier Excel Meteostat",
    filetypes=[("Fichiers Excel", "*.xlsx *.xls")]
)

# On demande le second fichier
fichier_excel_infoclimat = filedialog.askopenfilename(
    title="Choisir le fichier Excel Infoclimat",
    filetypes=[("Fichiers Excel", "*.xlsx *.xls")]
)

# Correction de la condition : s'assurer que les deux fichiers ont été sélectionnés
if not fichier_excel_meteostat or not fichier_excel_infoclimat:
    raise FileNotFoundError("Au moins un fichier n'a pas été sélectionné. Veuillez sélectionner les deux fichiers.")

# --- 2. Lecture et Préparation des données ---

### Lecture et traitement du fichier METEOSTAT (Data 1)
# Lis les colonnes nécessaires Meteostat : A = Date, B = Température, F = vent, H = pression
# Remarque : La colonne A (Date) est souvent lue par défaut comme un type string ou datetime.
try:
    data_meteostat = pd.read_excel(
        fichier_excel_meteostat, 
        usecols="A,B,F,H", 
        skiprows=1 # On saute l'en-tête (ligne 1)
        # nrows a été retiré pour lire toutes les lignes disponibles, à ajuster si nécessaire
    )
    # Renommer les colonnes pour une utilisation claire et cohérente
    data_meteostat.columns = ["Date", "Température moyenne", "Vent moyen", "Pression moyenne"]
    
    # Conversion de la colonne Date en format datetime
    # On utilise format='mixed' pour une meilleure robustesse
    data_meteostat["Date"] = pd.to_datetime(data_meteostat["Date"], errors="coerce", format='mixed')
    
    # Suppression des lignes avec Date manquante/invalide
    data_meteostat = data_meteostat.dropna(subset=["Date"])
    
    # Conversion des colonnes numériques (pour gérer les types de données non numériques qui pourraient causer des erreurs de tracé)
    for col in ["Température moyenne", "Vent moyen", "Pression moyenne"]:
        # 'coerce' va transformer les valeurs non convertibles en NaN
        data_meteostat[col] = pd.to_numeric(data_meteostat[col], errors='coerce')
        
    # La ligne suivante n'est plus nécessaire si l'Infoclimat ne contient pas d'humidité
    # On peut ajouter une colonne NaN pour l'humidité si on en a besoin pour une structure commune
    # data_meteostat["Humidité"] = np.nan 

except Exception as e:
    print(f"Erreur lors de la lecture ou du traitement du fichier Meteostat : {e}")
    # On lève une erreur ou on quitte si la lecture échoue
    raise

### Lecture et traitement du fichier INFOCLIMAT (Data 2)
# Lis les colonnes nécessaires Infoclimat : D = Date, F = Température, G = Pression, I = Humidité, L = Vent
# Remarque : Les indices D, F, G, I, L correspondent aux colonnes 4, 6, 7, 9, 12.
# [...]
### Lecture et traitement du fichier INFOCLIMAT (Data 2)
# Lis les colonnes nécessaires Infoclimat : D = Date, F = Température, G = Pression, I = Humidité, L = Vent
try:
    data_infoclimat = pd.read_excel(
        fichier_excel_infoclimat, 
        usecols="C,F,G,I,L", 
        skiprows=6 # On saute les premières lignes de l'en-tête (ligne 1 à 6)
    )
    
    # Renommer les colonnes
    data_infoclimat.columns = ["Date", "Température moyenne", "Pression moyenne", "Humidité", "Vent moyen"]
    
    # Conversion de la colonne Date en format datetime
    data_infoclimat["Date"] = pd.to_datetime(data_infoclimat["Date"], errors="coerce", format='mixed')
    
    # Suppression des lignes avec Date manquante/invalide
    data_infoclimat = data_infoclimat.dropna(subset=["Date"])

    # Conversion des colonnes numériques
    for col in ["Température moyenne", "Pression moyenne", "Humidité", "Vent moyen"]:
        data_infoclimat[col] = pd.to_numeric(data_infoclimat[col], errors='coerce')
        
    # =========================================================================
    # 🌟 NOUVEAU : Calcul de la moyenne journalière
    # =========================================================================
    
    # 1. Créer une colonne pour la date SANS l'heure
    data_infoclimat['Date_Jour'] = data_infoclimat['Date'].dt.normalize()
    
    # 2. Grouper par la date du jour et calculer la moyenne pour chaque colonne
    data_infoclimat_journalier = data_infoclimat.groupby('Date_Jour').agg({
        'Température moyenne': 'mean',
        'Pression moyenne': 'mean',
        'Humidité': 'mean',
        'Vent moyen': 'mean'
    }).reset_index()
    
    # 3. Remplacer le DataFrame original par le DataFrame moyenné
    # Ceci garantit que les tracés utiliseront les données journalières moyennées
    data_infoclimat = data_infoclimat_journalier.rename(columns={'Date_Jour': 'Date'})
    
    # =========================================================================
        
except Exception as e:
    print(f"Erreur lors de la lecture ou du traitement du fichier Infoclimat : {e}")
    raise
# [...]
    
    # Renommer les colonnes (attention à l'ordre D,F,G,I,L dans le fichier)
    data_infoclimat.columns = ["Date", "Température moyenne", "Pression moyenne", "Humidité", "Vent moyen"]
    
    # Conversion de la colonne Date en format datetime
    # La **GROSSE ERREUR** était ici : data2["Date"] = pd.to_datetime(data["Date"], errors="coerce")
    # qui utilisait la colonne "Date" du premier DataFrame (data) au lieu de data2.
    data_infoclimat["Date"] = pd.to_datetime(data_infoclimat["Date"], errors="coerce", format='mixed')
    
    # Suppression des lignes avec Date manquante/invalide
    data_infoclimat = data_infoclimat.dropna(subset=["Date"])

    # Conversion des colonnes numériques
    for col in ["Température moyenne", "Pression moyenne", "Humidité", "Vent moyen"]:
        data_infoclimat[col] = pd.to_numeric(data_infoclimat[col], errors='coerce')
        
except Exception as e:
    print(f"Erreur lors de la lecture ou du traitement du fichier Infoclimat : {e}")
    raise

# S'assurer que les dataframes ne sont pas vides après le nettoyage
if data_meteostat.empty or data_infoclimat.empty:
    raise ValueError("Un des DataFrames est vide après le nettoyage des dates.")

# --- 3. Tracé et Superposition des données ---
# S'assurer que les DataFrames sont triés par date pour un tracé correct
data_meteostat = data_meteostat.sort_values(by="Date")
data_infoclimat = data_infoclimat.sort_values(by="Date")

plt.figure(figsize=(15, 10)) # Taille de la figure ajustée pour une meilleure lisibilité

# Subplot 1 : Température
plt.subplot(2, 2, 1)
# Tracé des deux sources sur le même axe
plt.plot(data_meteostat["Date"], data_meteostat["Température moyenne"], label="Meteostat", color='tab:red', linestyle='-')
plt.plot(data_infoclimat["Date"], data_infoclimat["Température moyenne"], label="Infoclimat", color='darkred', linestyle='--')
plt.title("Température moyenne")
plt.xlabel("Date")
plt.ylabel("Température (°C)")
plt.legend() # Ajouter la légende pour identifier les sources
plt.grid(True)

# Subplot 2 : Pression Moyenne
plt.subplot(2, 2, 2)
# Tracé des deux sources sur le même axe
plt.plot(data_meteostat["Date"], data_meteostat["Pression moyenne"], label="Meteostat", color='tab:purple', linestyle='-')
plt.plot(data_infoclimat["Date"], data_infoclimat["Pression moyenne"], label="Infoclimat", color='indigo', linestyle='--')
plt.title("Pression moyenne")
plt.xlabel("Date")
plt.ylabel("Pression (hPa)")
plt.legend()
plt.grid(True)

# Subplot 3 : Humidité
plt.subplot(2, 2, 3)
# Attention : Vous n'aviez pas de colonne "Humidité" dans data_meteostat (votre premier fichier)
# Si Meteostat ne fournit pas d'humidité, seul le tracé Infoclimat sera pertinent, ou Meteostat tracera des NaN.
plt.plot(data_infoclimat["Date"], data_infoclimat["Humidité"], label="Infoclimat", color='tab:blue', linestyle='-')
# Si data_meteostat avait l'humidité :
# plt.plot(data_meteostat["Date"], data_meteostat["Humidité"], label="Meteostat", color='darkblue', linestyle='--')
plt.title("Humidité relative")
plt.xlabel("Date")
plt.ylabel("Humidité (%)")
plt.legend()
plt.grid(True)

# Subplot 4 : Vent moyen
plt.subplot(2, 2, 4)
# Tracé des deux sources sur le même axe
plt.plot(data_meteostat["Date"], data_meteostat["Vent moyen"], label="Meteostat", color='tab:green', linestyle='-')
plt.plot(data_infoclimat["Date"], data_infoclimat["Vent moyen"], label="Infoclimat", color='darkgreen', linestyle='--')
plt.title("Vent moyen")
plt.xlabel("Date")
plt.ylabel("Vent (m/s)")
plt.legend()
plt.grid(True)

plt.tight_layout() # Ajuster automatiquement les sous-graphiques pour qu'ils s'ajustent à la figure
plt.show()