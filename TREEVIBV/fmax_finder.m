%% ============================================================
%  analyse_cluster_script.m
%  Analyse multi-pics + suivi du cluster modal au fil du temps
% ============================================================
clearvars; close all; clc;
%% PARAMÈTRES UTILISATEUR
mat_folder = "H:\Home\Documents\aaaTREEVIB\downsampling_filtered";
fs_target = 25;                 
f_search = [0.05 0.6];           
Np = 10;                      % nombre max de pics significatifs
window_size=20;

t_start = 50;                    
t_end = 250;                    

min_window = 512;
max_window = 2048;
%% ==============================================
% LISTER LES FICHIERS
% ==============================================
files = dir(fullfile(mat_folder, "*.mat"));
if isempty(files)
    error("Aucun .mat trouvé dans %s", mat_folder);
end
fprintf("Fichiers trouvés : %d\n", numel(files));
%% ==============================================
%  STOCKAGE RÉSULTATS
% ==============================================
all_dates      = NaT(numel(files),1);
all_fcenter    = nan(numel(files),1);
all_fwidth     = nan(numel(files),1);
all_centroid   = nan(numel(files),1);   % <=== AJOUTE POUR LA CENTROÏDE
%% ==============================================
%  BOUCLE PRINCIPALE
% ==============================================
for k = 1:numel(files)
    fname = fullfile(files(k).folder, files(k).name);
    fprintf("[%d/%d] %s\n", k, numel(files), files(k).name);
    %% Chargement du fichier
    try
        S = load(fname);
    catch
        warning("Impossible de charger %s", files(k).name);
        continue;
    end
    %% Date
    dt = NaT;
    if isfield(S,'date_heure')
        try
            dt = datetime(S.date_heure);
        catch
            dt = NaT;
        end
    end
    if isnat(dt)
        dt = parse_filename_date(files(k).name);
    end
    all_dates(k) = dt;
    %% Extraction canal
    if isfield(S, 'Arbre1_Nord')
        x = S.Arbre1_Nord(:);
    else
        warning("Aucun canal 'Arbre1_Nord' dans %s", files(k).name);
        continue;
    end
    %% Fréquence d’échantillonnage
    if isfield(S,'fs')
        fs0 = double(S.fs);
    else
        fs0 = 1000;
        warning("fs manquant dans %s → fs=1000 utilisé", files(k).name);
    end
    %% RESAMPLING
    x = resample(x, fs_target, fs0);
    
    %% TRONCATURE TEMPORELLE (AJOUT)
    % Calcul des indices
    idx_start = round(t_start * fs_target) + 1; 
    idx_end = round(t_end * fs_target);
    % Assurer que les indices sont valides
    Nx_resampled = length(x);
    if idx_start >= idx_end || idx_start > Nx_resampled || idx_end < 1
        warning("Plage de temps demandée (%.1f-%.1f s) invalide ou dépasse le signal de longueur %.1f s dans %s", ...
                t_start, t_end, Nx_resampled/fs_target, files(k).name);
        continue; % Passe au fichier suivant si la plage est invalide
    end
    idx_start = max(1, idx_start);
    idx_end = min(Nx_resampled, idx_end);
    % Troncature du signal
    x = x(idx_start:idx_end);
    
    %% Prétraitement
    x = detrend(x,'linear');
    [b,a] = butter(4, 1/(fs_target/2), 'low');
    x = filtfilt(b,a,x);

    %% PSD par pwelch robuste
    Nx = length(x);
    Nw = min(max_window, max(min_window, floor(Nx/2)));
    win  = hamming(Nw);
    ov   = floor(Nw/2);
    nfft = 2^nextpow2(Nw);
    [Pxx,f] = pwelch(x, win, ov, nfft, fs_target);

    %% Sélection bande cluster
    m  = f >= f_search(1) & f <= f_search(2);
    ff = f(m);
    PP = Pxx(m);
    if numel(ff) < 5
        warning("Plage trop courte dans %s", files(k).name);
        continue;
    end
    %% ============================
    %   DETECTION MULTIPICS
    %% ============================
    [pks,locs,w,prom] = findpeaks(PP, ff, ...
        'MinPeakProminence', max(PP)*0.02, ...
        'MinPeakDistance', 0.02, ...
        'SortStr','descend');
    if isempty(pks)
        warning("Aucun pic détecté dans %s", files(k).name);
        continue;
    end
    Nkeep = min(Np, numel(pks));
    fpeaks = locs(1:Nkeep);
    amp    = pks(1:Nkeep);
    % tri croissant
    [fpeaks,is] = sort(fpeaks);
    amp = amp(is);
    %% ============================
    %   CLUSTER : CENTRE + LARGEUR
    %% ============================

    f_center = sum(fpeaks .* amp) / sum(amp);
    f_width  = fpeaks(end) - fpeaks(1);
    all_fcenter(k) = f_center;
    all_fwidth(k)  = f_width;
    %% ============================
    %   CENTROIDE SPECTRALE
    %% ============================

    centroid = sum(ff .* PP) / sum(PP);
    all_centroid(k) = centroid;
    fprintf(" → centre = %.3f Hz | width = %.3f Hz | centroid = %.3f Hz\n", ...
            f_center, f_width, centroid);
end

%% Lissage pour les tracés
valid = ~isnat(all_dates);

% all_fcenter_smoothed=smoothdata(all_fcenter(valid),window_size);
% all_centroid_smoothed=smoothdata(all_centroid(valid),window_size);

%% =======================================================
%  TRACÉS
% =======================================================

%% fréquence centrale
figure; hold on; grid on;
plot(all_dates(valid), all_fcenter(valid), 'x', 'LineWidth', 1.5);hold on
% plot(all_dates(valid), all_fcenter_smoothed,'r-')
ylabel("Fréquence centrale [Hz]");
xlabel("Date");
title("Centre du cluster");
%% largeur du cluster
figure; hold on; grid on;
plot(all_dates(valid), all_fwidth(valid), 'x', 'LineWidth', 1.5)
ylabel("Largeur [Hz]");
xlabel("Date");
title("Largeur du cluster modal");
%% centroïde spectrale
figure; hold on; grid on;
plot(all_dates(valid), all_centroid(valid), 'x', 'LineWidth', 1.5);hold on
% plot(all_dates(valid), all_centroid_smoothed,'r-')
ylabel("Centroïde spectrale [Hz]");
xlabel("Date");
title("Centroïde spectrale dans la bande");
disp("=== FIN ===");
%% ============================================
% FONCTION : parse date du nom
%% ============================================
function dt = parse_filename_date(fname)
    dt = NaT;
    [~,name,~] = fileparts(fname);
    tok = regexp(name, '(\d{4})(\d{2})(\d{2})[_-]?(\d{2})(\d{2})(\d{2})', 'tokens');
    if ~isempty(tok)
        v = tok{1};
        dt = datetime(str2double(v{1}), str2double(v{2}), str2double(v{3}), ...
                      str2double(v{4}), str2double(v{5}), str2double(v{6}));
    end
end