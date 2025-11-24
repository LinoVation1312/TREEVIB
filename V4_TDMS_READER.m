% Lecture et traitement de TOUS les fichiers TDMS d'un dossier avec filtrage
clear all; close all; clc;

%% ========== PARAMÈTRES CONFIGURABLES ==========
% Paramètres du dossier
tdms_folder = "V:\MONITORING_ARBRES\PEUPLIERS a partir du 03 mars (a mettre a jour chaque mois)";

% Paramètres du signal
fs_original = 1000; 
fs = 25;             

% Paramètres de filtrage
f_dc = 0.1;         % HP [Hz]
fc = 0.4;            % LP [Hz]
hp_taps = 2001;      % taps FIR HP
lp_order = 2;        % Ordre du filtre Butterworth

% Paramètres de pwelch
pwelch_window = 300*fs;
pwelch_overlap = 0 * pwelch_window;
pwelch_nfft = 2^nextpow2(pwelch_window);
pwelch_freq_range = [0 1]; 

% Organisation par arbre
trees = struct(...
    'Arbre1', {{'Arbre1_Nord', 'Arbre1_Est'}}, ...
    'Arbre2', {{'Arbre2_Nord', 'Arbre2_Est'}} ...
);


pwelch_freq_range = [0 1];
dynamic_range_dB = 30; 

% ================================================

%% Vérification du dossier
if ~isfolder(tdms_folder)
    error("Le dossier %s est introuvable.", tdms_folder);
end

%% Scan des fichiers TDMS
fprintf('\n=== Scan des fichiers TDMS dans le dossier ===\n');
tdms_files = dir(fullfile(tdms_folder, '*.tdms'));

if isempty(tdms_files)
    error('Aucun fichier TDMS trouvé dans le dossier %s', tdms_folder);
end

fprintf('Nombre de fichiers TDMS trouvés : %d\n', length(tdms_files));

% Trier par nom (ordre chronologique si nommage cohérent)
[~, sort_idx] = sort({tdms_files.name});
tdms_files = tdms_files(sort_idx);

% Afficher la liste
for i = 1:length(tdms_files)
    fprintf('  [%d] %s (%.2f MB)\n', i, tdms_files(i).name, ...
            tdms_files(i).bytes/1024/1024);
end

%% Analyse des groupes disponibles dans tous les fichiers
fprintf('\n=== Analyse des groupes TDMS ===\n');
all_groups = {};
valid_files = struct('name', {}, 'path', {}, 'groups', {}, 'channels', {});

for i = 1:length(tdms_files)
    file_path = fullfile(tdms_files(i).folder, tdms_files(i).name);
    
    try
        info = tdmsinfo(file_path);
        
        if ~isprop(info, 'ChannelList')
            fprintf('  ATTENTION : %s - pas de ChannelList, ignoré\n', tdms_files(i).name);
            continue;
        end
        
        chanlist = info.ChannelList;
        if ~all(ismember({'ChannelGroupName','ChannelName'}, chanlist.Properties.VariableNames))
            fprintf('  ATTENTION : %s - format ChannelList invalide, ignoré\n', tdms_files(i).name);
            continue;
        end
        
        file_groups = unique(chanlist.ChannelGroupName);
        all_groups = [all_groups; file_groups];
        
        % Stocker les infos du fichier valide
        valid_files(end+1).name = tdms_files(i).name;
        valid_files(end).path = file_path;
        valid_files(end).groups = file_groups;
        valid_files(end).channels = chanlist;
        
        fprintf('  OK : %s - %d groupe(s)\n', tdms_files(i).name, length(file_groups));
        
    catch ME
        fprintf('  ERREUR : %s - %s\n', tdms_files(i).name, ME.message);
    end
end

if isempty(valid_files)
    error('Aucun fichier TDMS valide trouvé.');
end

fprintf('\nFichiers valides : %d/%d\n', length(valid_files), length(tdms_files));

%% Consolidation des groupes
unique_groups = unique(all_groups);
fprintf('\n=== Groupes TDMS disponibles (tous fichiers confondus) ===\n');
disp(unique_groups);

%% Sélection du groupe
try
    [idx, ok] = listdlg('PromptString','Sélectionner un groupe:', ...
                        'SelectionMode','single', ...
                        'ListString', cellstr(unique_groups), ...
                        'Name','Sélection du groupe');
    if ~ok, error('Sélection annulée.'); end
    selected_group = unique_groups(idx);
catch
    fprintf('Interface indisponible. Sélection via console.\n');
    for gi = 1:numel(unique_groups)
        fprintf('  [%d] %s\n', gi, unique_groups{gi});
    end
    gi = input('Entrez l''indice du groupe à analyser: ');
    if ~(isnumeric(gi) && gi>=1 && gi<=numel(unique_groups))
        error('Indice invalide.');
    end
    selected_group = unique_groups(gi);
end

fprintf('\nGroupe sélectionné: %s\n', selected_group);

%% Identification des canaux du groupe (intersection de tous les fichiers)
fprintf('\n=== Identification des canaux communs ===\n');
common_channels = {};

for i = 1:length(valid_files)
    if ~ismember(selected_group, valid_files(i).groups)
        continue;
    end
    
    chanlist = valid_files(i).channels;
    in_group = strcmp(chanlist.ChannelGroupName, selected_group);
    file_channels = chanlist.ChannelName(in_group);
    
    if isempty(common_channels)
        common_channels = file_channels;
    else
        common_channels = intersect(common_channels, file_channels, 'stable');
    end
end

if isempty(common_channels)
    error('Aucun canal commun trouvé pour le groupe %s', selected_group);
end

fprintf('Canaux communs dans le groupe %s :\n', selected_group);
disp(common_channels);

%% Lecture et concaténation des données de tous les fichiers
fprintf('\n=== Lecture et concaténation des fichiers ===\n');
concatenated_data = struct();

for ch_idx = 1:length(common_channels)
    channel_name = common_channels{ch_idx};
    concatenated_data.(channel_name) = [];
end

for i = 1:length(valid_files)
    if ~ismember(selected_group, valid_files(i).groups)
        fprintf('  Skip : %s (groupe absent)\n', valid_files(i).name);
        continue;
    end
    
    fprintf('  Lecture : %s ... ', valid_files(i).name);
    
    try
        data_table = tdmsread(valid_files(i).path, ...
                              ChannelGroupName = selected_group, ...
                              ChannelNames = common_channels);
        
        % Déballer si nécessaire
        if iscell(data_table) && numel(data_table) == 1 && istable(data_table{1})
            data_table = data_table{1};
        end
        
        if ~istable(data_table)
            fprintf('ERREUR (format inattendu)\n');
            continue;
        end
        
        % Concaténer chaque canal
        for ch_idx = 1:length(common_channels)
            channel_name = common_channels{ch_idx};
            if ismember(channel_name, data_table.Properties.VariableNames)
                channel_data = double(data_table.(channel_name)(:));
                concatenated_data.(channel_name) = [concatenated_data.(channel_name); channel_data];
            end
        end
        
        fprintf('OK (%d points)\n', height(data_table));
        
    catch ME
        fprintf('ERREUR (%s)\n', ME.message);
    end
end

fprintf('\n=== Statistiques des données concaténées ===\n');
for ch_idx = 1:length(common_channels)
    channel_name = common_channels{ch_idx};
    n_points = length(concatenated_data.(channel_name));
    duration = n_points / fs_original;
    fprintf('  %s : %d points (%.2f s = %.2f min)\n', ...
            channel_name, n_points, duration, duration/60);
end

%% Séparer canaux temporels et DSP
is_dsp = endsWith(common_channels, "_DSP");
channels_time = common_channels(~is_dsp);
channels_dsp  = common_channels(is_dsp);

if isempty(channels_time)
    error('Aucun canal temporel présent dans le groupe.');
end

%% Affichage des paramètres
fprintf('\n=== Paramètres de traitement ===\n');
fprintf('Fréquence originale : %d Hz\n', fs_original);
fprintf('Fréquence après resampling : %d Hz\n', fs);
fprintf('Taille de fenêtre pwelch : %d\n', pwelch_window);
fprintf('Recouvrement : %d échantillons\n', pwelch_overlap);
fprintf('NFFT : %d\n', pwelch_nfft);
fprintf('Résolution fréquentielle : %.4f Hz\n', fs/pwelch_nfft);

%% ========== TRAITEMENT ET STOCKAGE AVEC NOMS EXPLICITES ==========
DATA = struct();

fprintf('\n=== Traitement de tous les canaux ===\n');

for ch_idx = 1:length(channels_time)
    channel_name = channels_time{ch_idx};
    
    if ~isfield(concatenated_data, channel_name)
        fprintf('ATTENTION : Canal %s non trouvé, ignoré.\n', channel_name);
        continue;
    end
    
    fprintf('  -> %s\n', channel_name);
    
    % Extraction des données concaténées à fs_original
    data_1000Hz = concatenated_data.(channel_name);
    N_original = length(data_1000Hz);
    
    % Resampling à fs Hz
    N_target = round(N_original * fs / fs_original);
    data = resample(data_1000Hz, N_target, N_original);
    clear data_1000Hz;
    
    fprintf('     Points à %d Hz : %d (durée: %.2f min)\n', ...
            fs, length(data), length(data) / fs / 60);
    
    N = length(data);
    t = (0:N-1) / fs;
    
    [Pxx, f_psd] = pwelch(data, pwelch_window, pwelch_overlap, pwelch_nfft, fs);
    
    data_bp = filtfilt_cascade(data, fs, f_dc, fc, hp_taps, lp_order);
    
    [Pxx_bp, f_psd_bp] = pwelch(data_bp, pwelch_window, pwelch_overlap, pwelch_nfft, fs);
    
    var_name = [lower(strrep(channel_name, ' ', '_')), '_time'];
    var_name = matlab.lang.makeValidName(var_name);
    
 
    DATA.(var_name) = struct(...
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

fprintf('\n=== Variables disponibles dans DATA ===\n');
available_vars = fieldnames(DATA);
for i = 1:length(available_vars)
    orig_name = DATA.(available_vars{i}).channel_name;
    fprintf('  DATA.%s  (canal: %s)\n', available_vars{i}, orig_name);
end

%% ========== AFFICHAGE DES FIGURES (Hors Fonction - Comparaison Arbres 1 et 2) ==========


colors_dir = [0 0.447 0.741; 0.85 0.33 0.1]; % Bleu pour Nord, Rouge pour Est

% Vérification rapide de la présence des données nécessaires pour éviter les erreurs
required_vars = {'arbre1_nord_time', 'arbre1_est_time', 'arbre2_nord_time', 'arbre2_est_time'};
for i = 1:length(required_vars)
    if ~isfield(DATA, required_vars{i})
        error('Données manquantes : DATA.%s n''est pas disponible. Vérifiez le nom des canaux dans la structure "trees".', required_vars{i});
    end
end

% Récupération des données pour une écriture plus concise
D1_N = DATA.arbre1_nord_time;
D1_E = DATA.arbre1_est_time;
D2_N = DATA.arbre2_nord_time;
D2_E = DATA.arbre2_est_time;

%% === FIGURE 1 : Brut (non filtré) ===
fprintf('\n=== Création de la Figure 1 : Brut (Temporel et PSD) ===\n');

figure('Name', 'Figure 1: Signaux Bruts et PSD Brutes', 'Position', [100, 100, 1400, 900]);
sgtitle(sprintf('Signaux Bruts - Groupe %s', selected_group), ...
        'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');

% --- Subplot 1: Signal temporel brut - Arbre 1 ---
subplot(2,2,1); hold on; grid on;
plot(D1_N.t, D1_N.data, 'Color', colors_dir(1,:), 'DisplayName', 'Arbre 1 Nord');
plot(D1_E.t, D1_E.data, 'Color', colors_dir(2,:), 'DisplayName', 'Arbre 1 Est');
title('Arbre 1 - Temporel Brut'); xlabel('t [s]'); ylabel('Amplitude');
legend('show', 'Location', 'best', 'Interpreter','none');

% --- Subplot 2: Signal temporel brut - Arbre 2 ---
subplot(2,2,2); hold on; grid on;
plot(D2_N.t, D2_N.data, 'Color', colors_dir(1,:), 'DisplayName', 'Arbre 2 Nord');
plot(D2_E.t, D2_E.data, 'Color', colors_dir(2,:), 'DisplayName', 'Arbre 2 Est');
title('Arbre 2 - Temporel Brut'); xlabel('t [s]'); ylabel('Amplitude');
legend('show', 'Location', 'best', 'Interpreter','none');

% --- Calcul de la plage dynamique pour les PSD brutes (Arbre 1) ---
Pxx_1N_dB = 10*log10(D1_N.Pxx);
Pxx_1E_dB = 10*log10(D1_E.Pxx);
idx_range_1 = D1_N.f_psd >= pwelch_freq_range(1) & D1_N.f_psd <= pwelch_freq_range(2);

max_overall_1_raw = max([max(Pxx_1N_dB(idx_range_1)); max(Pxx_1E_dB(idx_range_1))]);
ylim_raw_1 = [max_overall_1_raw - dynamic_range_dB, max_overall_1_raw + 5];

% --- Subplot 3: PSD brute - Arbre 1 ---
subplot(2,2,3); hold on; grid on;
plot(D1_N.f_psd, Pxx_1N_dB, 'Color', colors_dir(1,:), 'DisplayName', 'Arbre 1 Nord');
plot(D1_E.f_psd, Pxx_1E_dB, 'Color', colors_dir(2,:), 'DisplayName', 'Arbre 1 Est');
title('Arbre 1 - PSD Brute');
xlabel('f [Hz]'); ylabel('PSD [dB]');
xlim(pwelch_freq_range);
ylim(ylim_raw_1); % Application de la plage dynamique
legend('show', 'Location', 'best', 'Interpreter','none');


% --- Calcul de la plage dynamique pour les PSD brutes (Arbre 2) ---
Pxx_2N_dB = 10*log10(D2_N.Pxx);
Pxx_2E_dB = 10*log10(D2_E.Pxx);
idx_range_2 = D2_N.f_psd >= pwelch_freq_range(1) & D2_N.f_psd <= pwelch_freq_range(2);

max_overall_2_raw = max([max(Pxx_2N_dB(idx_range_2)); max(Pxx_2E_dB(idx_range_2))]);
ylim_raw_2 = [max_overall_2_raw - dynamic_range_dB, max_overall_2_raw + 5];

% --- Subplot 4: PSD brute - Arbre 2 ---
subplot(2,2,4); hold on; grid on;
plot(D2_N.f_psd, Pxx_2N_dB, 'Color', colors_dir(1,:), 'DisplayName', 'Arbre 2 Nord');
plot(D2_E.f_psd, Pxx_2E_dB, 'Color', colors_dir(2,:), 'DisplayName', 'Arbre 2 Est');
title('Arbre 2 - PSD Brute');
xlabel('f [Hz]'); ylabel('PSD [dB]');
xlim(pwelch_freq_range);
ylim(ylim_raw_2); % Application de la plage dynamique
legend('show', 'Location', 'best', 'Interpreter','none');


%% === FIGURE 2 : Filtré ===
fprintf('\n=== Création de la Figure 2 : Filtré (Temporel et PSD) ===\n');

figure('Name', 'Figure 2: Signaux Filtrés et PSD Filtrées', 'Position', [150, 150, 1400, 900]);
sgtitle(sprintf('Signaux Filtrés (%s-%s Hz) - Groupe %s', num2str(f_dc), num2str(fc), selected_group), ...
        'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');

% --- Subplot 1: Signal temporel filtré - Arbre 1 ---
subplot(2,2,1); hold on; grid on;
plot(D1_N.t, D1_N.data_bp, 'Color', colors_dir(1,:), 'DisplayName', 'Arbre 1 Nord');
plot(D1_E.t, D1_E.data_bp, 'Color', colors_dir(2,:), 'DisplayName', 'Arbre 1 Est');
title('Arbre 1 - Temporel Filtré'); xlabel('t [s]'); ylabel('Amplitude');
legend('show', 'Location', 'best', 'Interpreter','none');

% --- Subplot 2: Signal temporel filtré - Arbre 2 ---
subplot(2,2,2); hold on; grid on;
plot(D2_N.t, D2_N.data_bp, 'Color', colors_dir(1,:), 'DisplayName', 'Arbre 2 Nord');
plot(D2_E.t, D2_E.data_bp, 'Color', colors_dir(2,:), 'DisplayName', 'Arbre 2 Est');
title('Arbre 2 - Temporel Filtré'); xlabel('t [s]'); ylabel('Amplitude');
legend('show', 'Location', 'best', 'Interpreter','none');

% --- Calcul de la plage dynamique pour les PSD filtrées (Arbre 1) ---
Pxx_bp_1N_dB = 10*log10(D1_N.Pxx_bp);
Pxx_bp_1E_dB = 10*log10(D1_E.Pxx_bp);

max_overall_1_bp = max([max(Pxx_bp_1N_dB(idx_range_1)); max(Pxx_bp_1E_dB(idx_range_1))]);
ylim_bp_1 = [max_overall_1_bp - dynamic_range_dB, max_overall_1_bp + 5];

% --- Subplot 3: PSD filtrée - Arbre 1 ---
subplot(2,2,3); hold on; grid on;
plot(D1_N.f_psd_bp, Pxx_bp_1N_dB, 'Color', colors_dir(1,:), 'DisplayName', 'Arbre 1 Nord');
plot(D1_E.f_psd_bp, Pxx_bp_1E_dB, 'Color', colors_dir(2,:), 'DisplayName', 'Arbre 1 Est');
title('Arbre 1 - PSD Filtrée');
xlabel('f [Hz]'); ylabel('PSD [dB]');
xlim(pwelch_freq_range);
ylim(ylim_bp_1); % Application de la plage dynamique
legend('show', 'Location', 'best', 'Interpreter','none');


% --- Calcul de la plage dynamique pour les PSD filtrées (Arbre 2) ---
Pxx_bp_2N_dB = 10*log10(D2_N.Pxx_bp);
Pxx_bp_2E_dB = 10*log10(D2_E.Pxx_bp);

max_overall_2_bp = max([max(Pxx_bp_2N_dB(idx_range_2)); max(Pxx_bp_2E_dB(idx_range_2))]);
ylim_bp_2 = [max_overall_2_bp - dynamic_range_dB, max_overall_2_bp + 5];

% --- Subplot 4: PSD filtrée - Arbre 2 ---
subplot(2,2,4); hold on; grid on;
plot(D2_N.f_psd_bp, Pxx_bp_2N_dB, 'Color', colors_dir(1,:), 'DisplayName', 'Arbre 2 Nord');
plot(D2_E.f_psd_bp, Pxx_bp_2E_dB, 'Color', colors_dir(2,:), 'DisplayName', 'Arbre 2 Est');
title('Arbre 2 - PSD Filtrée');
xlabel('f [Hz]'); ylabel('PSD [dB]');
xlim(pwelch_freq_range);
ylim(ylim_bp_2); % Application de la plage dynamique
legend('show', 'Location', 'best', 'Interpreter','none');



%% Tracés orbitaux
% Créer la figure
fig = figure('Position', [100 100 800 600]);

% Créer les subplots et stocker les handles des tracés
ax1 = subplot(1,2,1);
h1 = plot(ax1, D1_E.data_bp, D1_N.data_bp, 'b.');
xlabel('D1 E (m)');
ylabel('D1 N (m)');
title('Signal D1');
grid on;
axis tight;

ax2 = subplot(1,2,2);
h2 = plot(ax2, D2_E.data_bp, D2_N.data_bp, 'r.');
xlabel('D2 E (m)');
ylabel('D2 N (m)');
title('Signal D2');
grid on;
axis tight;

% Créer les sliders
% Slider pour t_min
uicontrol('Style', 'text', ...
    'Position', [50 50 50 20], ...
    'String', 't_min:');

slider_tmin = uicontrol('Style', 'slider', ...
    'Position', [110 50 300 20], ...
    'Min', 0, 'Max', duration, 'Value', 0, ...
    'SliderStep', [1/duration, 10/duration]);

txt_tmin = uicontrol('Style', 'text', ...
    'Position', [420 50 50 20], ...
    'String', '0.0 s');

uicontrol('Style', 'text', ...
    'Position', [50 20 50 20], ...
    'String', 't_max:');

slider_tmax = uicontrol('Style', 'slider', ...
    'Position', [110 20 300 20], ...
    'Min', 0.1, 'Max', duration, 'Value', duration, ...
    'SliderStep', [1/duration, 10/duration]);

txt_tmax = uicontrol('Style', 'text', ...
    'Position', [420 20 50 20], ...
    'String', sprintf('%.1f s', duration));

% Stocker les handles dans UserData de la figure
handles_data.slider_tmin = slider_tmin;
handles_data.slider_tmax = slider_tmax;
handles_data.txt_tmin = txt_tmin;
handles_data.txt_tmax = txt_tmax;
handles_data.h1 = h1;
handles_data.h2 = h2;
handles_data.ax1 = ax1;
handles_data.ax2 = ax2;
handles_data.D1_E = D1_E;
handles_data.D1_N = D1_N;
handles_data.D2_E = D2_E;
handles_data.D2_N = D2_N;
handles_data.fs = fs;
set(fig, 'UserData', handles_data);

% Assigner les callbacks APRÈS avoir stocké les données
set(slider_tmin, 'Callback', {@updatePlots, fig});
set(slider_tmax, 'Callback', {@updatePlots, fig});

% Initialiser les tracés
updatePlots([], [], fig);

%% ========== FIN DU SCRIPT PRINCIPAL ==========
% Mettez cette fonction à la fin de votre fichier V3_TDMS_READER.m

function updatePlots(~, ~, fig)
    % Récupérer les handles depuis UserData
    handles_data = get(fig, 'UserData');
    
    slider_tmin = handles_data.slider_tmin;
    slider_tmax = handles_data.slider_tmax;
    txt_tmin = handles_data.txt_tmin;
    txt_tmax = handles_data.txt_tmax;
    h1 = handles_data.h1;
    h2 = handles_data.h2;
    ax1 = handles_data.ax1;
    ax2 = handles_data.ax2;
    D1_E = handles_data.D1_E;
    D1_N = handles_data.D1_N;
    D2_E = handles_data.D2_E;
    D2_N = handles_data.D2_N;
    fs = handles_data.fs;
    
    % Récupérer les valeurs des sliders
    t_min = get(slider_tmin, 'Value');
    t_max = get(slider_tmax, 'Value');
    
    % S'assurer que t_min < t_max
    if t_min >= t_max
        t_min = t_max - 0.1;
        set(slider_tmin, 'Value', t_min);
    end
    
    % Mettre à jour les textes
    set(txt_tmin, 'String', sprintf('%.1f s', t_min));
    set(txt_tmax, 'String', sprintf('%.1f s', t_max));
    
    % Calculer les indices correspondants
    idx_min = round(t_min * fs) + 1;
    idx_max = round(t_max * fs) + 1;
    
    % Limiter les indices aux dimensions des données
    idx_min = max(1, idx_min);
    idx_max = min(length(D1_E.data_bp), idx_max);
    
    % Mettre à jour les données des tracés
    set(h1, 'XData', D1_E.data_bp(idx_min:idx_max), ...
            'YData', D1_N.data_bp(idx_min:idx_max));
    set(h2, 'XData', D2_E.data_bp(idx_min:idx_max), ...
            'YData', D2_N.data_bp(idx_min:idx_max));
    
    % Réajuster les axes
    axis(ax1, 'tight');
    axis(ax2, 'tight');
end
%% ========= Fonctions auxiliaires =========
function y = filtfilt_cascade(x, fs, f_dc, fc, hp_taps, lp_order)
    % Filtre passe-haut FIR (soustraction DC) + passe-bas Butterworth, appliqués avec filtfilt
    
    % HP FIR (type I, atténuation DC)
    nyq = fs/2;
    if f_dc <= 0
        y_hp = x;
    else
        Wn_hp = f_dc/nyq;
        b_hp = fir1(hp_taps-1, Wn_hp, 'high', hamming(hp_taps), 'scale');
        y_hp = filtfilt(b_hp, 1, x);
    end

    % LP Butterworth
    if fc >= nyq
        y = y_hp;
        return;
    end
    Wn_lp = fc/nyq;
    [b_lp, a_lp] = butter(lp_order, Wn_lp, 'low');
    y = filtfilt(b_lp, a_lp, y_hp);
end


