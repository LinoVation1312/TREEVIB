from nptdms import TdmsFile
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.signal import butter, filtfilt, firwin,resample
from scipy.fft import fft, fftfreq
from scipy.io import wavfile

plt.close('all')

def bandpass_dc_lp(x, fs, f_hp=0.01, f_lp=0.5, hp_taps=10001, lp_order=2):

    nyq = 0.5 * fs

    # --- Passe-haut FIR pour DC / très basse fréquence ---
    taps = firwin(hp_taps, f_hp / nyq, pass_zero=False)
    y = filtfilt(taps, [1.0], x)

    # --- Passe-bas Butterworth ---
    b, a = butter(lp_order, f_lp / nyq, btype='low')
    y = filtfilt(b, a, y)

    return y

tdms_path = "H:\\Home\\Telechargements\\mesure_2025_05_16_15_36_15.tdms"


if not os.path.exists(tdms_path):
    raise FileNotFoundError(f"Le fichier {tdms_path} est introuvable.")


tdms_file = TdmsFile.read(tdms_path)


print("\n=== Informations TDMS ===")
print(f"Groupes trouvés : {tdms_file.groups()}")
print("--------------------------")


for group in tdms_file.groups():
    print(f"\nGROUPE : {group.name}")
    for channel in group.channels():
        print(f" Canal : {channel.name} — {len(channel)} points")


group_name = tdms_file.groups()[0].name # premier groupe
channel_name = tdms_file.groups()[0].channels()[0].name # premier canal
data = tdms_file[group_name][channel_name].data


df = tdms_file.as_dataframe(time_index=False)
print("\nAperçu des données :")
print(df.head())



N=len(data);
fs=1000;
t=np.arange(N)/fs;

fc = 0.5
f_dc=0.05


data_bp = bandpass_dc_lp(data, fs=fs, f_hp=f_dc, f_lp=fc)



Y = fft(data)
f = fftfreq(N, 1/fs)  # axe fréquence en Hz

plt.figure(figsize=(10,6))
plt.subplot(2,2,1)
plt.plot(t, data)
plt.title(f"{group_name} - {channel_name} (signal brut)")
plt.xlabel("Time [s]")
plt.ylabel("Acc. [m/s²]")
plt.grid(True)

plt.subplot(2,2,2)
plt.plot(f[:N//2], 2.0/N * np.abs(Y[:N//2]))  # FFT positive seulement
plt.title("FFT (amplitude)")
plt.xlabel("Frequency [Hz]")
plt.ylabel("Amplitude")
plt.grid(True)
plt.xlim(0, 1)
plt.tight_layout()
plt.show()
plt.xlim(0, 1)

# --- FFT du signal filtré ---
Yf = fft(data_bp)

plt.subplot(2,2,3)
plt.plot(t, data_bp)
plt.title(f"{group_name} - {channel_name} (signal filtré)")
plt.xlabel("Time [s]")
plt.ylabel("Acc. [m/s²]")
plt.grid(True)

plt.subplot(2,2,4)
plt.plot(f[:N//2], 2.0/N * np.abs(Yf[:N//2]))
plt.title("FFT (filtré)")
plt.xlabel("Frequency [Hz]")
plt.ylabel("Amplitude")
plt.grid(True)
plt.xlim(0, 1)
plt.tight_layout()
plt.show()



# === DOWNSAMPLING à 20 Hz ===
fs_target = 20  # Hz


N_target = int(N * fs_target / fs)
data_20hz = resample(data, N_target)
data_bp_20hz = resample(data_bp, N_target)



data_20hz_wav = data_20hz
data_bp_20hz_wav = data_bp_20hz

# === ENREGISTREMENT EN WAV ===
output_dir = "H:\\Home\\Telechargements"
wav_raw = os.path.join(output_dir, "signal_brut_20Hz.wav")
wav_filtered = os.path.join(output_dir, "signal_filtre_20Hz.wav")

wavfile.write(wav_raw, fs_target, data_20hz_wav)
wavfile.write(wav_filtered, fs_target, data_bp_20hz_wav)

print(f"\n=== Fichiers WAV créés ===")
print(f"Signal brut (20 Hz) : {wav_raw}")
print(f"Signal filtré (20 Hz) : {wav_filtered}")
print(f"Durée : {len(data_20hz) / fs_target:.2f} s")
print(f"Nombre d'échantillons : {len(data_20hz)}")
