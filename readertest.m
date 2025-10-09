% Analyse de 2 fichiers WAV avec affichage temporel et FFT
clear; close all; clc;

%% Chargement des fichiers WAV
[file1, path1] = uigetfile('*.wav', 'Sélectionner le premier fichier WAV');
if isequal(file1, 0)
    disp('Aucun fichier sélectionné');
    return;
end

[file2, path2] = uigetfile('*.wav', 'Sélectionner le deuxième fichier WAV');
if isequal(file2, 0)
    disp('Aucun fichier sélectionné');
    return;
end

% Lecture des fichiers
[data1, fs1] = audioread(fullfile(path1, file1));
[data2, fs2] = audioread(fullfile(path2, file2));

% Si stéréo, prendre seulement le premier canal
if size(data1, 2) > 1
    data1 = data1(:, 1);
end
if size(data2, 2) > 1
    data2 = data2(:, 1);
end

%% Création des vecteurs temps
N1 = length(data1);
N2 = length(data2);
t1 = (0:N1-1) / fs1;
t2 = (0:N2-1) / fs2;

%% Calcul des FFT
Y1 = fft(data1);
Y2 = fft(data2);

f1 = (0:N1-1) * (fs1/N1);
f2 = (0:N2-1) * (fs2/N2);

% Spectre amplitude (unilatéral)
P1 = abs(Y1/N1);
P1 = P1(1:N1/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f1_half = f1(1:N1/2+1);

P2 = abs(Y2/N2);
P2 = P2(1:N2/2+1);
P2(2:end-1) = 2*P2(2:end-1);
f2_half = f2(1:N2/2+1);

%% Affichage
figure('Position', [100, 100, 1400, 800]);

% Signal 1 - Temporel
subplot(2, 2, 1);
plot(t1, data1, 'b', 'LineWidth', 1);
grid on;
title(['Signal temporel - ', file1], 'Interpreter', 'none');
xlabel('Temps [s]');
ylabel('Amplitude');
xlim([0 max(t1)]);

% Signal 1 - FFT
subplot(2, 2, 2);
plot(f1_half, P1, 'b', 'LineWidth', 1);
grid on;
title(['FFT - ', file1], 'Interpreter', 'none');
xlabel('Fréquence [Hz]');
ylabel('Amplitude');
xlim([0 1]);

% Signal 2 - Temporel
subplot(2, 2, 3);
plot(t2, data2, 'r', 'LineWidth', 1);
grid on;
title(['Signal temporel - ', file2], 'Interpreter', 'none');
xlabel('Temps [s]');
ylabel('Amplitude');
xlim([0 max(t2)]);

% Signal 2 - FFT
subplot(2, 2, 4);
plot(f2_half, P2, 'r', 'LineWidth', 1);
grid on;
title(['FFT - ', file2], 'Interpreter', 'none');
xlabel('Fréquence [Hz]');
ylabel('Amplitude');
xlim([0 1]);

%% Affichage des informations
fprintf('\n=== Informations Fichier 1 ===\n');
fprintf('Nom: %s\n', file1);
fprintf('Fréquence échantillonnage: %d Hz\n', fs1);
fprintf('Nombre échantillons: %d\n', N1);
fprintf('Durée: %.2f s\n', N1/fs1);

fprintf('\n=== Informations Fichier 2 ===\n');
fprintf('Nom: %s\n', file2);
fprintf('Fréquence échantillonnage: %d Hz\n', fs2);
fprintf('Nombre échantillons: %d\n', N2);
fprintf('Durée: %.2f s\n', N2/fs2);