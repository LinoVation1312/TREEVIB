% Lecture et traitement de TOUS les fichiers TDMS d'un dossier avec filtrage
clear all; close all; clc;

%% ========== PARAMÈTRES CONFIGURABLES ==========
% Paramètres du dossier
tdms_folder = "V:\MONITORING_ARBRES\PEUPLIERS a partir du 03 mars (a mettre a jour chaque mois)";

% Paramètres du signal
fs_original = 1000; 
fs = 40;             

% Paramètres de filtrage
f_dc = 0.1;         % HP [Hz]
fc = 0.5;            % LP [Hz]
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

%% Adapter la structure trees aux canaux disponibles
tree_names = fieldnames(trees);
for ti = 1:numel(tree_names)
    chs = trees.(tree_names{ti});
    trees.(tree_names{ti}) = intersect(chs, channels_time, 'stable');
end

% Supprimer les arbres vides
tree_names = fieldnames(trees);
empty_mask = cellfun(@(t) isempty(trees.(t)), tree_names);
trees = rmfield(trees, tree_names(empty_mask));
tree_names = fieldnames(trees);

if isempty(tree_names)
    trees = struct('Groupe', cellstr(channels_time(:)'));
    tree_names = fieldnames(trees);
    fprintf('Aucun arbre prédéfini valide. Analyse générique par canaux.\n');
end

colors = [0 0.447 0.741; 0.85 0.33 0.1];

%% Traitement par arbre
output_root = fullfile(pwd);

for tree_idx = 1:length(tree_names)
    tree_name = tree_names{tree_idx};
    channels = trees.(tree_name);

    fprintf('\n========================================\n');
    fprintf('Traitement de : %s (groupe %s)\n', tree_name, selected_group);
    fprintf('========================================\n');

    data_processed = struct();

    % Traitement de chaque canal
    for dir_idx = 1:length(channels)
        channel_name = channels{dir_idx};

        if ~isfield(concatenated_data, channel_name)
            fprintf('ATTENTION : Canal %s non trouvé, ignoré.\n', channel_name);
            continue;
        end

        fprintf('\n  -> %s\n', channel_name);

        % Extraction des données concaténées à fs_original
        data_1000Hz = concatenated_data.(channel_name);
        N_original = length(data_1000Hz);

        % Resampling à fs Hz
        N_target = round(N_original * fs / fs_original);
        data = resample(data_1000Hz, N_target, N_original);
        clear data_1000Hz;

        fprintf('     Points à %d Hz : %d\n', fs, length(data));
        fprintf('     Durée : %.2f s (%.2f min)\n', length(data) / fs, length(data) / fs / 60);

        % Calcul du temps
        N = length(data);
        t = (0:N-1) / fs;

        % PSD brute
        [Pxx, f_psd] = pwelch(data, pwelch_window, pwelch_overlap, pwelch_nfft, fs);

        % Filtrage
        data_bp = filtfilt_cascade(data, fs, f_dc, fc, hp_taps, lp_order);

        % PSD filtrée
        [Pxx_bp, f_psd_bp] = pwelch(data_bp, pwelch_window, pwelch_overlap, pwelch_nfft, fs);

        % Stockage
        data_processed.(channel_name) = struct(...
            't', t, ...
            'data', data, ...
            'data_bp', data_bp, ...
            'Pxx', Pxx, ...
            'f_psd', f_psd, ...
            'Pxx_bp', Pxx_bp, ...
            'f_psd_bp', f_psd_bp ...
        );
    end

    %% Figure
    figure('Position', [100, 100, 1400, 900]);
    figure_title = sprintf('Analyse %s (groupe %s) - %d fichiers - fs = %d Hz', ...
                           tree_name, selected_group, length(valid_files), fs);
    sgtitle(figure_title, 'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');

    % Tracé 1: signal brut
    subplot(2,2,1); hold on; grid on;
    for dir_idx = 1:length(channels)
        channel_name = channels{dir_idx};
        if ~isfield(data_processed, channel_name), continue; end
        d = data_processed.(channel_name);
        c = colors(1 + mod(dir_idx-1, size(colors,1)), :);
        plot(d.t, d.data, 'Color', c, 'DisplayName', channel_name);
    end
    xlabel('Temps [s]'); ylabel('Amplitude'); title('Signal brut'); legend('Interpreter','none');

    % Tracé 2: signal filtré
    subplot(2,2,2); hold on; grid on;
    for dir_idx = 1:length(channels)
        channel_name = channels{dir_idx};
        if ~isfield(data_processed, channel_name), continue; end
        d = data_processed.(channel_name);
        c = colors(1 + mod(dir_idx-1, size(colors,1)), :);
        plot(d.t, d.data_bp, 'Color', c, 'DisplayName', channel_name);
    end
    xlabel('Temps [s]'); ylabel('Amplitude'); title('Signal filtré (HP+LP)'); legend('Interpreter','none');

    % Tracé 3: PSD brute
    subplot(2,2,3); hold on; grid on;
    for dir_idx = 1:length(channels)
        channel_name = channels{dir_idx};
        if ~isfield(data_processed, channel_name), continue; end
        d = data_processed.(channel_name);
        c = colors(1 + mod(dir_idx-1, size(colors,1)), :);
        plot(d.f_psd, d.Pxx, 'Color', c, 'DisplayName', channel_name);
        xlim(pwelch_freq_range);
    end
    xlabel('Fréquence [Hz]'); ylabel('PSD [m/s²/Hz]'); title('PSD brute'); legend('Interpreter','none');

    % Tracé 4: PSD filtrée
    subplot(2,2,4); hold on; grid on;
    for dir_idx = 1:length(channels)
        channel_name = channels{dir_idx};
        if ~isfield(data_processed, channel_name), continue; end
        d = data_processed.(channel_name);
        c = colors(1 + mod(dir_idx-1, size(colors,1)), :);
        plot(d.f_psd_bp, d.Pxx_bp, 'Color', c, 'DisplayName', channel_name);
        xlim(pwelch_freq_range);
    end
    xlabel('Fréquence [Hz]'); ylabel('[m/s²/Hz]'); title('PSD filtrée'); legend('Interpreter','none');
end

fprintf('\n=== Terminé ===\n');

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