# -*- coding: utf-8 -*-
"""
Script : tracer plusieurs variables météo depuis un fichier Excel
"""

import pandas as pd
import matplotlib.pyplot as plt
from tkinter import Tk, filedialog

# --- 1. Sélection du fichier Excel ---
root = Tk()
root.withdraw()  # cache la fenêtre principale Tkinter
fichier_excel = filedialog.askopenfilename(
    title="Choisir le fichier Excel",
    filetypes=[("Fichiers Excel", "*.xlsx *.xls")]
)

if not fichier_excel:
    raise FileNotFoundError("Aucun fichier sélectionné.")

# --- 2. Lecture des données ---
# Lis toutes les colonnes nécessaires : D = Date, F = Température, G = Pression, H = Humidité, I = Vent
data = pd.read_excel(fichier_excel, usecols="A,B,F,H", skiprows=1, nrows=214)

# Renommer les colonnes pour plus de clarté
data.columns = ["Date", "Température moyenne", "Vent moyen", "Pression moyenne"]

# Conversion de la colonne Date en format datetime
data["Date"] = pd.to_datetime(data["Date"], errors="coerce")

# Suppression des lignes avec données manquantes
data = data.dropna(subset=["Date"])

# --- 3. Tracé ---
plt.figure(figsize=(12,8))

# # Subplot 1 : toutes les variables sur le même graphique
# plt.subplot(2,2,1)
# plt.plot(data["Date"], data["Température moyenne"], label="Température (°C)")
# plt.plot(data["Date"], data["Pression moyenne"], label="Pression (hPa)")
# plt.plot(data["Date"], data["Humidité"], label="Humidité (%)")
# plt.plot(data["Date"], data["Vent moyen"], label="Vent (m/s)")
# plt.title("Évolution des paramètres météo")
# plt.xlabel("Date")
# plt.ylabel("Valeurs mesurées")
# plt.legend()
# plt.grid(True)

# Subplot 1 : température seule
plt.subplot(2,2,1)
plt.plot(data["Date"], data["Température moyenne"], color='tab:red')
plt.title("Température moyenne")
plt.xlabel("Date")
plt.ylabel("Température (°C)")
plt.grid(True)

# Subplot 2 : Pression Moyenne seule
plt.subplot(2,2,2)
plt.plot(data["Date"], data["Pression moyenne"], color='tab:purple')
plt.title("Pression moyenne")
plt.xlabel("Date")
plt.ylabel("Pression (hPa)")
plt.grid(True)

# # Subplot 3 : humidité
# plt.subplot(2,2,3)
# plt.plot(data["Date"], data["Humidité"], color='tab:blue')
# plt.title("Humidité relative")
# plt.xlabel("Date")
# plt.ylabel("Humidité (%)")
# plt.grid(True)

# Subplot 4 : vent moyen
plt.subplot(2,2,4)
plt.plot(data["Date"], data["Vent moyen"], color='tab:green')
plt.title("Vent moyen")
plt.xlabel("Date")
plt.ylabel("Vent (m/s)")
plt.grid(True)

plt.tight_layout()
plt.show()
