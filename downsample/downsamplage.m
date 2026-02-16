% Extraction et Downsampling des signaux temporels (TDMS vers MAT)
clear all; close all; clc;

%% ========== PARAMÈTRES CONFIGURABLES ==========
% Paramètres du dossier
tdms_folder = "V:\MONITORING_ARBRES\PEUPLIERS a partir du 03 mars (a mettre a jour chaque mois)";

% Paramètres du signal
fs_original = 1000; 
fs_target = 5;             

% Organisation par arbre (utilisé pour identifier les canaux)
trees = struct(...
    'Arbre1', {{'Arbre1_Nord', 'Arbre1_Est'}}, ...
    'Arbre2', {{'Arbre2_Nord', 'Arbre2_Est'}} ...
);

% Les paramètres de filtrage, pwelch, fmax, etc. ont été supprimés.
% ================================================

%% 1. Vérification du dossier source et Sélection du Dossier de Sortie via GUI 📁

if ~isfolder(tdms_folder)
    error("Le dossier TDMS %s est introuvable.", tdms_folder);
end

fprintf('\n=== Sélection du dossier de sortie ===\n');
output_folder_temp = uigetdir(tdms_folder, 'Sélectionnez le dossier pour enregistrer les signaux downsamplés (5 Hz)');

if output_folder_temp == 0 % L'utilisateur a annulé
    error('Opération annulée par l''utilisateur.');
end

output_folder = output_folder_temp;
fprintf('Dossier de sortie sélectionné : %s\n', output_folder);

%% 2. Scan des fichiers TDMS et Préparation des groupes
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

%% Consolidation des groupes et Extraction des dates
unique_groups = unique(all_groups);
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


%% 3. Traitement, Downsampling et Sauvegarde du Signal ⬇️

fprintf('\n=== Traitement des groupes (Extraction et Downsampling à %d Hz) ===\n', fs_target);

for group_idx = 1:length(unique_groups_sorted)
    selected_group = unique_groups_sorted{group_idx};
    current_date = group_dates_sorted(group_idx);
    
    if isnat(current_date)
        fprintf('[%d/%d] Skip : %s (date non parsable)\n', group_idx, length(unique_groups_sorted), selected_group);
        continue;
    end
    
    date_str = datestr(current_date, 'yyyymmdd_HHMMSS');
    fprintf('[%d/%d] Traitement : %s (%s)\n', group_idx, length(unique_groups_sorted), ...
            selected_group, datestr(current_date, 'dd/mm/yyyy HH:MM:SS'));
    
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
    
    %% Downsampling et Sauvegarde pour chaque canal
    
    % Extraire les noms des arbres depuis la structure 'trees'
    tree_names = fieldnames(trees);
    
    for ti = 1:length(tree_names)
        tree_name = tree_names{ti};
        directions = trees.(tree_name);
        
        for di = 1:length(directions)
            channel_name = directions{di};
            
            if ~isfield(concatenated_data, channel_name) || isempty(concatenated_data.(channel_name))
                continue;
            end
            
            % Rééchantillonnage (Downsampling)
            data_original_fs = concatenated_data.(channel_name);
            N_original = length(data_original_fs);
            % Calculer le rapport pour l'interpolation/décimation
            p = fs_target;
            q = fs_original;
            
            data_downsampled = resample(data_original_fs, p, q);
            
            % Création du nom de fichier
            % Format: yyyymmdd_HHMMSS_ChannelName.mat
            file_name = sprintf('%s_%s.mat', date_str, channel_name);
            save_path = fullfile(output_folder, file_name);
            
            % Sauvegarde des données
            signal = data_downsampled;
            date_heure = current_date;
            fs = fs_target;
            save(save_path, 'signal', 'date_heure', 'fs');
            
            fprintf('  -> %s : %d pts (original %d) sauvegardés dans %s\n', ...
                    channel_name, length(data_downsampled), N_original, file_name);
        end
    end
end

fprintf('\n=== Extraction et Downsampling terminés. Données sauvegardées dans le dossier choisi. ===\n');


%% 4. Fonctions auxiliaires

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