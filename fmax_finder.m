% Analyse temporelle de fmax - TDMS multiples fichiers avec filtrage
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

% Paramètres de recherche de fmax
fmax_search_range = [0.05 0.8];  % Plage de recherche du pic [Hz]

% Paramètres de lissage temporel
smooth_method = 'movmean';  % 'movmean', 'movmedian', 'gaussian', 'lowess', 'rlowess'
smooth_window = 7;          % Taille de fenêtre pour le lissage (en jours)

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

% Trier par nom (ordre chronologique)
[~, sort_idx] = sort({tdms_files.name});
tdms_files = tdms_files(sort_idx);

%% Analyse des groupes disponibles dans tous les fichiers
fprintf('\n=== Analyse des groupes TDMS ===\n');
all_groups = {};
valid_files = struct('name', {}, 'path', {}, 'groups', {}, 'channels', {});

for i = 1:length(tdms_files)
    file_path = fullfile(tdms_files(i).folder, tdms_files(i).name);
    
    try
        info = tdmsinfo(file_path);
        
        if ~isprop(info, 'ChannelList')
            continue;
        end
        
        chanlist = info.ChannelList;
        if ~all(ismember({'ChannelGroupName','ChannelName'}, chanlist.Properties.VariableNames))
            continue;
        end
        
        file_groups = unique(chanlist.ChannelGroupName);
        all_groups = [all_groups; file_groups];
        
        valid_files(end+1).name = tdms_files(i).name;
        valid_files(end).path = file_path;
        valid_files(end).groups = file_groups;
        valid_files(end).channels = chanlist;
        
    catch ME
        fprintf('  ERREUR : %s - %s\n', tdms_files(i).name, ME.message);
    end
end

if isempty(valid_files)
    error('Aucun fichier TDMS valide trouvé.');
end

fprintf('Fichiers valides : %d/%d\n', length(valid_files), length(tdms_files));

%% Consolidation des groupes
unique_groups = unique(all_groups);
fprintf('\n=== Groupes TDMS disponibles ===\n');
fprintf('Total : %d groupes\n', length(unique_groups));

%% Extraction des dates depuis les noms de groupes
group_dates = NaT(size(unique_groups));
for i = 1:length(unique_groups)
    group_name = unique_groups{i};
    date_parsed = parse_group_date(group_name);
    if ~isnat(date_parsed)
        group_dates(i) = date_parsed;
    end
end

% Trier par date
[group_dates_sorted, sort_idx] = sort(group_dates);
unique_groups_sorted = unique_groups(sort_idx);

fprintf('Plage temporelle : %s -> %s\n', ...
        datestr(min(group_dates(~isnat(group_dates)))), ...
        datestr(max(group_dates(~isnat(group_dates)))));

%% Initialisation de la structure de résultats
% Structure: results(tree)(direction).dates / .fmax / .group_names
tree_names = fieldnames(trees);
results = struct();

for ti = 1:length(tree_names)
    tree_name = tree_names{ti};
    directions = trees.(tree_name);
    
    for di = 1:length(directions)
        dir_name = directions{di};
        results.(tree_name).(dir_name) = struct(...
            'dates', NaT(0,1), ...
            'fmax', [], ...
            'group_names', {{}} ...
        );
    end
end

%% Traitement de chaque groupe
fprintf('\n=== Traitement des groupes ===\n');

for group_idx = 1:length(unique_groups_sorted)
    selected_group = unique_groups_sorted{group_idx};
    current_date = group_dates_sorted(group_idx);
    
    if isnat(current_date)
        fprintf('[%d/%d] Skip : %s (date non parsable)\n', group_idx, length(unique_groups_sorted), selected_group);
        continue;
    end
    
    fprintf('[%d/%d] Traitement : %s (%s)\n', group_idx, length(unique_groups_sorted), ...
            selected_group, datestr(current_date, 'dd/mm/yyyy HH:MM'));
    
    %% Identification des canaux communs
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
        fprintf('  -> Aucun canal commun, ignoré\n');
        continue;
    end
    
    %% Lecture et concaténation
    concatenated_data = struct();
    for ch_idx = 1:length(common_channels)
        channel_name = common_channels{ch_idx};
        concatenated_data.(channel_name) = [];
    end
    
    for i = 1:length(valid_files)
        if ~ismember(selected_group, valid_files(i).groups)
            continue;
        end
        
        try
            data_table = tdmsread(valid_files(i).path, ...
                                  ChannelGroupName = selected_group, ...
                                  ChannelNames = common_channels);
            
            if iscell(data_table) && numel(data_table) == 1 && istable(data_table{1})
                data_table = data_table{1};
            end
            
            if ~istable(data_table)
                continue;
            end
            
            for ch_idx = 1:length(common_channels)
                channel_name = common_channels{ch_idx};
                if ismember(channel_name, data_table.Properties.VariableNames)
                    channel_data = double(data_table.(channel_name)(:));
                    concatenated_data.(channel_name) = [concatenated_data.(channel_name); channel_data];
                end
            end
            
        catch ME
            continue;
        end
    end
    
    %% Traitement par arbre/direction
    for ti = 1:length(tree_names)
        tree_name = tree_names{ti};
        directions = trees.(tree_name);
        
        for di = 1:length(directions)
            channel_name = directions{di};
            
            if ~isfield(concatenated_data, channel_name) || isempty(concatenated_data.(channel_name))
                continue;
            end
            
            % Resampling
            data_1000Hz = concatenated_data.(channel_name);
            N_original = length(data_1000Hz);
            N_target = round(N_original * fs / fs_original);
            data = resample(data_1000Hz, N_target, N_original);
            
            % Vérification de la longueur minimale
            if length(data) < pwelch_window
                fprintf('  -> %s : signal trop court (%d pts < %d), ignoré\n', ...
                        channel_name, length(data), pwelch_window);
                continue;
            end
            
            % Filtrage
            data_bp = filtfilt_cascade(data, fs, f_dc, fc, hp_taps, lp_order);
            
            % PSD filtrée
            [Pxx_bp, f_psd_bp] = pwelch(data_bp, pwelch_window, pwelch_overlap, pwelch_nfft, fs);
            
            % Recherche de fmax dans la plage définie
            freq_mask = (f_psd_bp >= fmax_search_range(1)) & (f_psd_bp <= fmax_search_range(2));
            [~, max_idx] = max(Pxx_bp(freq_mask));
            f_in_range = f_psd_bp(freq_mask);
            fmax = f_in_range(max_idx);
            
            % Stockage
            results.(tree_name).(channel_name).dates(end+1) = current_date;
            results.(tree_name).(channel_name).fmax(end+1) = fmax;
            results.(tree_name).(channel_name).group_names{end+1} = selected_group;
            
            fprintf('  -> %s : fmax = %.4f Hz\n', channel_name, fmax);
        end
    end
end

%% Moyennage par jour
fprintf('\n=== Moyennage par jour ===\n');

results_daily = struct();

for ti = 1:length(tree_names)
    tree_name = tree_names{ti};
    directions = trees.(tree_name);
    
    for di = 1:length(directions)
        dir_name = directions{di};
        
        if ~isfield(results.(tree_name), dir_name) || isempty(results.(tree_name).(dir_name).dates)
            results_daily.(tree_name).(dir_name) = struct('dates', NaT(0,1), 'fmax', []);
            continue;
        end
        
        dates_raw = results.(tree_name).(dir_name).dates;
        fmax_raw = results.(tree_name).(dir_name).fmax;
        
        % Conversion en jours (sans heure)
        dates_day = dateshift(dates_raw, 'start', 'day');
        unique_days = unique(dates_day);
        
        fmax_daily = zeros(size(unique_days));
        
        for day_idx = 1:length(unique_days)
            day_mask = (dates_day == unique_days(day_idx));
            fmax_daily(day_idx) = mean(fmax_raw(day_mask));
        end
        
        results_daily.(tree_name).(dir_name).dates = unique_days;
        results_daily.(tree_name).(dir_name).fmax = fmax_daily;
        
        fprintf('%s - %s : %d mesures -> %d jours\n', tree_name, dir_name, ...
                length(dates_raw), length(unique_days));
    end
end

%% Tracé de l'évolution temporelle
fprintf('\n=== Génération des graphiques ===\n');

colors_dir = [0 0.447 0.741; 0.85 0.33 0.1];  % Bleu, Orange

figure('Position', [100, 100, 1600, 900]);
sgtitle('Évolution temporelle de fmax (signal filtré)', 'FontSize', 16, 'FontWeight', 'bold');

plot_idx = 1;

for ti = 1:length(tree_names)
    tree_name = tree_names{ti};
    directions = trees.(tree_name);
    
    for di = 1:length(directions)
        dir_name = directions{di};
        
        subplot(2, 2, plot_idx);
        hold on; grid on;
        
        if ~isfield(results_daily.(tree_name), dir_name) || isempty(results_daily.(tree_name).(dir_name).dates)
            title(sprintf('%s - %s (pas de données)', tree_name, dir_name), 'Interpreter', 'none');
            plot_idx = plot_idx + 1;
            continue;
        end
        
        dates = results_daily.(tree_name).(dir_name).dates;
        fmax_vals = results_daily.(tree_name).(dir_name).fmax;
        
        % Tracer les points bruts
        color = colors_dir(1 + mod(di-1, size(colors_dir, 1)), :);
        plot(dates, fmax_vals, 'o', 'Color', color, 'MarkerSize', 6, ...
             'DisplayName', 'Mesures journalières', 'MarkerFaceColor', color);
        
        % Lissage
        if length(fmax_vals) > smooth_window
            fmax_smooth = smoothdata(fmax_vals, smooth_method, smooth_window);
            plot(dates, fmax_smooth, '-', 'Color', color, 'LineWidth', 2, ...
                 'DisplayName', sprintf('Lissage (%s, %d jours)', smooth_method, smooth_window));
        end
        
        xlabel('Date');
        ylabel('fmax [Hz]');
        title(sprintf('%s - %s', tree_name, dir_name), 'Interpreter', 'none');
        legend('Location', 'best');
        datetick('x', 'dd/mm', 'keepticks', 'keeplimits');
        xtickangle(45);
        
        plot_idx = plot_idx + 1;
    end
end

fprintf('\n=== Terminé ===\n');

%% ========= Fonctions auxiliaires =========

function date_parsed = parse_group_date(group_name)
    % Parse les noms de groupe au format : 2025_07_23_00_00_07
    % Retourne datetime ou NaT si échec
    
    date_parsed = NaT;
    
    try
        % Extraction par expression régulière
        tokens = regexp(group_name, '(\d{4})_(\d{2})_(\d{2})_(\d{2})_(\d{2})_(\d{2})', 'tokens');
        
        if ~isempty(tokens) && length(tokens{1}) == 6
            year = str2double(tokens{1}{1});
            month = str2double(tokens{1}{2});
            day = str2double(tokens{1}{3});
            hour = str2double(tokens{1}{4});
            minute = str2double(tokens{1}{5});
            second = str2double(tokens{1}{6});
            
            date_parsed = datetime(year, month, day, hour, minute, second);
        end
    catch
        % Retourne NaT en cas d'erreur
    end
end

function y = filtfilt_cascade(x, fs, f_dc, fc, hp_taps, lp_order)
    % Filtre passe-haut FIR (soustraction DC) + passe-bas Butterworth
    
    nyq = fs/2;
    
    % HP FIR
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