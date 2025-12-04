% Lecture et traitement des fichiers MAT temporels avec filtrage
clear all; close all; clc;

%% ========== PARAMÈTRES CONFIGURABLES ==========
mat_folder = "H:\Home\Documents\aaaTREEVIB\downsampling_filtered"; % REMPLACEZ PAR VOTRE DOSSIER
% mat_folder = "V:\MONITORING_ARBRES\PEUPLIERS a partir du 03 mars (a mettre a jour chaque mois)";

% Paramètres du signal
fs_target = 25;             % Fréquence d'échantillonnage après resampling [Hz]
default_fs = 1000;          % Fréquence par défaut si 'fs' n'est pas trouvé dans le .mat

% Paramètres de filtrage (pour les tracés temporels/orbitaux)
f_dc = 0.1;                 % Coupure HP [Hz]
fc = 0.4;                   % Coupure LP [Hz]
hp_taps = 2001;             % Taps FIR HP
lp_order = 2;               % Ordre du filtre Butterworth

% Paramètres de pwelch
% Ces valeurs doivent être cohérentes avec fs_target (ex: 300*25 = 7500 points)
pwelch_window = 300 * fs_target;
pwelch_overlap = 0;
pwelch_nfft = 2^nextpow2(pwelch_window);
pwelch_freq_range = [0 1]; 
dynamic_range_dB = 30; 

% Noms des canaux à lire (doivent exister dans les fichiers .mat)
% Utilisation de la même structure que l'original pour la compatibilité
trees = struct(...
    'Arbre1', {{'Arbre1_Nord', 'Arbre1_Est'}}, ...
    'Arbre2', {{'Arbre2_Nord', 'Arbre2_Est'}} ...
);
channel_names = unique([trees.Arbre1, trees.Arbre2]); % {'Arbre1_Nord', 'Arbre1_Est', 'Arbre2_Nord', 'Arbre2_Est'}

% Couleurs pour les tracés Nord/Est
colors_dir = [0 0.447 0.741; 0.85 0.33 0.1]; % Bleu pour Nord, Rouge pour Est

%% ================================================
%% Vérification du dossier
if ~isfolder(mat_folder)
    error("Le dossier %s est introuvable.", mat_folder);
end

%% Scan et Sélection des fichiers MAT
fprintf('\n=== Scan des fichiers MAT dans le dossier ===\n');
all_files = dir(fullfile(mat_folder, '*.mat'));
if isempty(all_files)
    error('Aucun fichier MAT trouvé dans le dossier %s', mat_folder);
end

% Trier par nom (ordre chronologique si nommage cohérent)
[~, sort_idx] = sort({all_files.name});
all_files = all_files(sort_idx);

% Préparer la liste pour la sélection
file_names_list = {all_files.name}';

% --- Sélection de(s) fichier(s) ---
fprintf('Nombre de fichiers MAT trouvés : %d\n', length(all_files));
try
    [selection_indices, ok] = listdlg('PromptString','Sélectionner les fichiers MAT à analyser:', ...
                                      'SelectionMode','multiple', ...
                                      'ListString', file_names_list, ...
                                      'Name','Sélection des Données');
    
    if ~ok || isempty(selection_indices)
        error('Sélection annulée ou aucun fichier sélectionné.');
    end
    
    mat_files = all_files(selection_indices);
    fprintf('Nombre de fichiers sélectionnés : %d\n', length(mat_files));

catch
    % Option de repli si l'interface listdlg échoue (ex: environnement sans affichage)
    warning('Interface graphique indisponible. Tous les fichiers seront traités.');
    mat_files = all_files;
end


%% ========== BOUCLE PRINCIPALE DE TRAITEMENT ==========
DATA_by_file = {}; % Pour stocker les résultats de chaque fichier

for k = 1:numel(mat_files)
    fname = fullfile(mat_files(k).folder, mat_files(k).name);
    file_name_only = mat_files(k).name;
    fprintf('\n[%d/%d] Traitement de %s\n', k, numel(mat_files), file_name_only);

    % --- Chargement du fichier ---
    try
        S = load(fname);
    catch ME
        warning('Impossible de charger %s : %s. Ignoré.', file_name_only, ME.message);
        continue;
    end

    % --- Date ---
    dt = parse_filename_date(file_name_only, S);
    if isnat(dt)
        warning('Date non trouvée pour %s. Ignoré.', file_name_only);
        continue;
    end

    % --- Fréquence d'échantillonnage ---
    fs0 = default_fs;
    if isfield(S,'fs')
        fs0 = double(S.fs);
    elseif isfield(S,'FS')
        fs0 = double(S.FS);
    end
    fprintf('  → Fs original : %d Hz\n', fs0);
    
    % Vérification rapide des canaux
    found_channels = intersect(channel_names, fieldnames(S), 'stable');
    if isempty(found_channels)
        warning('  → Aucun canal requis (%s) trouvé dans %s. Ignoré.', strjoin(channel_names, ', '), file_name_only);
        continue;
    end
    
    current_data = struct('date', dt, 'fs', fs_target, 'file_name', file_name_only);
    
    % --- Traitement canal par canal ---
    for ch_idx = 1:length(found_channels)
        channel_name = found_channels{ch_idx};
        
        % Extraction des données
        data_1000Hz = S.(channel_name)(:);
        N_original = length(data_1000Hz);
        
        % Resampling
        N_target = round(N_original * fs_target / fs0);
        data = resample(data_1000Hz, N_target, N_original);
        clear data_1000Hz;
        
        N = length(data);
        t = (0:N-1) / fs_target;
        
        fprintf('  → %s : %d points (%.2f min) après resampling.\n', ...
                channel_name, N, N / fs_target / 60);

        % PSD Brute
        [Pxx, f_psd] = pwelch(data, pwelch_window, pwelch_overlap, pwelch_nfft, fs_target);
        
        % Filtrage Passe-Bande
        data_bp = filtfilt_cascade(data, fs_target, f_dc, fc, hp_taps, lp_order);
        
        % PSD Filtrée
        [Pxx_bp, f_psd_bp] = pwelch(data_bp, pwelch_window, pwelch_overlap, pwelch_nfft, fs_target);
        
        % Stockage dans une structure unique par souci de lisibilité dans les figures
        var_name = [lower(strrep(channel_name, ' ', '_')), '_time'];
        var_name = matlab.lang.makeValidName(var_name);
        
        current_data.(var_name) = struct(...
            'channel_name', channel_name, ...
            't', t, ...
            'data', data, ...
            'data_bp', data_bp, ...
            'Pxx', Pxx, ...
            'f_psd', f_psd, ...
            'Pxx_bp', Pxx_bp, ...
            'f_psd_bp', f_psd_bp ...
        );
    end
    
    % Stockage des données de ce fichier
    DATA_by_file{k} = current_data;
    
    % --- Génération des figures pour ce fichier ---
    % Nécessite les 4 canaux pour les figures (adapter si moins de 4 canaux sont trouvés)
    required_vars = {'arbre1_nord_time', 'arbre1_est_time', 'arbre2_nord_time', 'arbre2_est_time'};
    has_all_channels = true;
    for rv = 1:length(required_vars)
        if ~isfield(current_data, required_vars{rv})
            has_all_channels = false;
            break;
        end
    end
    
    if has_all_channels
        fprintf('  → Génération des figures.\n');
        
        % Récupération des données pour une écriture plus concise
        A1_N = current_data.arbre1_nord_time;
        A1_E = current_data.arbre1_est_time;
        A2_N = current_data.arbre2_nord_time;
        A2_E = current_data.arbre2_est_time;
        
        % Appel aux fonctions de tracé (définies plus bas)
        plot_raw_data(A1_N, A1_E, A2_N, A2_E, dt, file_name_only, colors_dir, ...
                      pwelch_freq_range, dynamic_range_dB);
        plot_filtered_data(A1_N, A1_E, A2_N, A2_E, dt, file_name_only, colors_dir, ...
                           pwelch_freq_range, dynamic_range_dB, f_dc, fc);
        plot_orbitals(A1_N, A1_E, A2_N, A2_E, dt, file_name_only);
    else
        fprintf('  → ATTENTION: Les 4 canaux requis pour les figures ne sont pas présents. Figures ignorées.\n');
    end
end

fprintf('\n=== Traitement terminé. %d fichiers traités. ===\n', length(DATA_by_file));

%% ========== FONCTIONS DE TRACÉ ET DE TRAITEMENT ==========

function plot_raw_data(A1_N, A1_E, A2_N, A2_E, dt, file_name, colors_dir, pwelch_freq_range, dynamic_range_dB)
    % === FIGURE 1 : Brut ===
    figure('Name', ['Figure 1: Brut - ' file_name], 'Position', [100, 100, 1400, 900]);
    sgtitle(sprintf('Signaux Bruts - %s\n(%s)', file_name, datestr(dt)), ...
            'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Signal temporel brut - Arbre 1
    subplot(2,2,1); hold on; grid on;
    plot(A1_N.t, A1_N.data, 'Color', colors_dir(1,:), 'DisplayName', 'Arbre 1 Nord');
    plot(A1_E.t, A1_E.data, 'Color', colors_dir(2,:), 'DisplayName', 'Arbre 1 Est');
    title('Arbre 1 - Temporel Brut'); xlabel('t [s]'); ylabel('Amplitude [m/s²]');
    legend('show', 'Location', 'best', 'Interpreter','none');
    
    % Signal temporel brut - Arbre 2
    subplot(2,2,2); hold on; grid on;
    plot(A2_N.t, A2_N.data, 'Color', colors_dir(1,:), 'DisplayName', 'Arbre 2 Nord');
    plot(A2_E.t, A2_E.data, 'Color', colors_dir(2,:), 'DisplayName', 'Arbre 2 Est');
    title('Arbre 2 - Temporel Brut'); xlabel('t [s]'); ylabel('Amplitude [m/s²]');
    legend('show', 'Location', 'best', 'Interpreter','none');
    
    % Calcul dynamique PSD Arbre 1
    [Pxx_1N_dB, Pxx_1E_dB, idx_range_1, ylim_raw_1] = calculate_psd_limits(A1_N.Pxx, A1_E.Pxx, A1_N.f_psd, pwelch_freq_range, dynamic_range_dB);

    % PSD brute Arbre 1
    subplot(2,2,3); hold on; grid on;
    plot(A1_N.f_psd, Pxx_1N_dB, 'Color', colors_dir(1,:), 'DisplayName', 'Arbre 1 Nord');
    plot(A1_E.f_psd, Pxx_1E_dB, 'Color', colors_dir(2,:), 'DisplayName', 'Arbre 1 Est');
    title('Arbre 1 - PSD Brute');
    xlabel('f [Hz]'); ylabel('PSD [dB]');
    xlim(pwelch_freq_range);
    ylim(ylim_raw_1);
    legend('show', 'Location', 'best', 'Interpreter','none');
    
    % Calcul dynamique PSD Arbre 2
    [Pxx_2N_dB, Pxx_2E_dB, ~, ylim_raw_2] = calculate_psd_limits(A2_N.Pxx, A2_E.Pxx, A2_N.f_psd, pwelch_freq_range, dynamic_range_dB);
    
    % PSD brute - Arbre 2 
    subplot(2,2,4); hold on; grid on;
    plot(A2_N.f_psd, Pxx_2N_dB, 'Color', colors_dir(1,:), 'DisplayName', 'Arbre 2 Nord');
    plot(A2_E.f_psd, Pxx_2E_dB, 'Color', colors_dir(2,:), 'DisplayName', 'Arbre 2 Est');
    title('Arbre 2 - PSD Brute');
    xlabel('f [Hz]'); ylabel('PSD [dB]');
    xlim(pwelch_freq_range);
    ylim(ylim_raw_2);
    legend('show', 'Location', 'best', 'Interpreter','none');
end

function plot_filtered_data(A1_N, A1_E, A2_N, A2_E, dt, file_name, colors_dir, pwelch_freq_range, dynamic_range_dB, f_dc, fc)
    % === FIGURE 2 : Filtré ===
    figure('Name', ['Figure 2: Filtré - ' file_name], 'Position', [150, 150, 1400, 900]);
    sgtitle(sprintf('Signaux Filtrés (%s-%s Hz) - %s\n(%s)', num2str(f_dc), num2str(fc), file_name, datestr(dt)), ...
            'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Signal temporel filtré - Arbre 1
    subplot(2,2,1); hold on; grid on;
    plot(A1_N.t, A1_N.data_bp, 'Color', colors_dir(1,:), 'DisplayName', 'Arbre 1 Nord');
    plot(A1_E.t, A1_E.data_bp, 'Color', colors_dir(2,:), 'DisplayName', 'Arbre 1 Est');
    title('Arbre 1 - Temporel Filtré'); xlabel('t [s]'); ylabel('Amplitude [m/s²]');
    legend('show', 'Location', 'best', 'Interpreter','none');
    
    % Signal temporel filtré - Arbre 2
    subplot(2,2,2); hold on; grid on;
    plot(A2_N.t, A2_N.data_bp, 'Color', colors_dir(1,:), 'DisplayName', 'Arbre 2 Nord');
    plot(A2_E.t, A2_E.data_bp, 'Color', colors_dir(2,:), 'DisplayName', 'Arbre 2 Est');
    title('Arbre 2 - Temporel Filtré'); xlabel('t [s]'); ylabel('Amplitude [m/s²]');
    legend('show', 'Location', 'best', 'Interpreter','none');
    
    % Calcul dynamique PSD filtrées Arbre 1
    [Pxx_bp_1N_dB, Pxx_bp_1E_dB, idx_range_1, ylim_bp_1] = calculate_psd_limits(A1_N.Pxx_bp, A1_E.Pxx_bp, A1_N.f_psd_bp, pwelch_freq_range, dynamic_range_dB);
    
    % PSD filtrée - Arbre 1 
    subplot(2,2,3); hold on; grid on;
    plot(A1_N.f_psd_bp, Pxx_bp_1N_dB, 'Color', colors_dir(1,:), 'DisplayName', 'Arbre 1 Nord');
    plot(A1_E.f_psd_bp, Pxx_bp_1E_dB, 'Color', colors_dir(2,:), 'DisplayName', 'Arbre 1 Est');
    title('Arbre 1 - PSD Filtrée');
    xlabel('f [Hz]'); ylabel('PSD [dB]');
    xlim(pwelch_freq_range);
    ylim(ylim_bp_1);
    legend('show', 'Location', 'best', 'Interpreter','none');
    
    % Calcul dynamique PSD filtrées Arbre 2
    [Pxx_bp_2N_dB, Pxx_bp_2E_dB, ~, ylim_bp_2] = calculate_psd_limits(A2_N.Pxx_bp, A2_E.Pxx_bp, A2_N.f_psd_bp, pwelch_freq_range, dynamic_range_dB);
    
    % PSD filtrée - Arbre 2
    subplot(2,2,4); hold on; grid on;
    plot(A2_N.f_psd_bp, Pxx_bp_2N_dB, 'Color', colors_dir(1,:), 'DisplayName', 'Arbre 2 Nord');
    plot(A2_E.f_psd_bp, Pxx_bp_2E_dB, 'Color', colors_dir(2,:), 'DisplayName', 'Arbre 2 Est');
    title('Arbre 2 - PSD Filtrée');
    xlabel('f [Hz]'); ylabel('PSD [dB]');
    xlim(pwelch_freq_range);
    ylim(ylim_bp_2);
    legend('show', 'Location', 'best', 'Interpreter','none');
end

function plot_orbitals(A1_N, A1_E, A2_N, A2_E, dt, file_name)
    % === Tracés orbitaux ===
    fig = figure('Name', ['Figure 3: Orbitaux - ' file_name], 'Position', [100 100 800 600]);
    sgtitle(sprintf('Tracés Orbitaux Filtrés (%s)', datestr(dt)), ...
            'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');
            
    % Arbre 1
    ax1 = subplot(1,2,1);
    plot(ax1, A1_E.data_bp, A1_N.data_bp, 'b.');
    xlabel('D1 E (m/s²) [Filtré]'); % Mettre à jour l'unité (assumant accélération)
    ylabel('D1 N (m/s²) [Filtré]');
    title('Trajectoire A1 (Nord vs Est)');
    grid on;
    axis tight equal;
    
    % Arbre 2
    ax2 = subplot(1,2,2);
    plot(ax2, A2_E.data_bp, A2_N.data_bp, 'r.');
    xlabel('D2 E (m/s²) [Filtré]');
    ylabel('D2 N (m/s²) [Filtré]');
    title('Trajectoire A2 (Nord vs Est)');
    grid on;
    axis tight equal;
    
    % Pas de slider ici car l'analyse est sur un seul fichier
end

function dt = parse_filename_date(fname, S)
    % Tente de lire la date à partir du champ 'date_heure' puis du nom du fichier
    dt = NaT;
    if isfield(S,'date_heure')
        try
            dt = datetime(S.date_heure);
        catch
            dt = NaT;
        end
    end
    if isnat(dt)
        [~,name,~] = fileparts(fname);
        % Regex pour YYYYMMDD[_-]HHMMSS
        tok = regexp(name, '(\d{4})(\d{2})(\d{2})[_-]?(\d{2})(\d{2})(\d{2})', 'tokens');
        if ~isempty(tok)
            v = tok{1};
            dt = datetime(str2double(v{1}), str2double(v{2}), str2double(v{3}), ...
                          str2double(v{4}), str2double(v{5}), str2double(v{6}));
        end
    end
end

function y = filtfilt_cascade(x, fs, f_dc, fc, hp_taps, lp_order)
    % Filtre passe-haut FIR (soustraction DC) + passe-bas Butterworth, appliqués avec filtfilt
    % Identique à la fonction de votre script original
    nyq = fs/2;
    if f_dc <= 0
        y_hp = x;
    else
        Wn_hp = f_dc/nyq;
        b_hp = fir1(hp_taps-1, Wn_hp, 'high', hamming(hp_taps), 'scale');
        y_hp = filtfilt(b_hp, 1, x);
    end
    if fc >= nyq
        y = y_hp;
        return;
    end
    Wn_lp = fc/nyq;
    [b_lp, a_lp] = butter(lp_order, Wn_lp, 'low');
    y = filtfilt(b_lp, a_lp, y_hp);
end

function [Pxx_N_dB, Pxx_E_dB, idx_range, ylim_range] = calculate_psd_limits(Pxx_N, Pxx_E, f_psd, pwelch_freq_range, dynamic_range_dB)
    % Calcule les limites dynamiques pour les tracés PSD en dB
    Pxx_N_dB = 10*log10(Pxx_N);
    Pxx_E_dB = 10*log10(Pxx_E);
    idx_range = f_psd >= pwelch_freq_range(1) & f_psd <= pwelch_freq_range(2);
    
    % S'assurer qu'il y a des points dans la plage pour le calcul du max
    if any(idx_range)
        max_overall_raw = max([max(Pxx_N_dB(idx_range)); max(Pxx_E_dB(idx_range))]);
        ylim_range = [max_overall_raw - dynamic_range_dB, max_overall_raw + 5];
    else
        % Valeurs par défaut si la plage est vide
        max_overall_raw = max([max(Pxx_N_dB); max(Pxx_E_dB)]);
        ylim_range = [max_overall_raw - dynamic_range_dB, max_overall_raw + 5];
    end
end