% Concaténation des fichiers MAT par date/heure
% Regroupe tous les canaux d'une même date/heure dans un seul fichier .mat
clear all; close all; clc;

%% ========== PARAMÈTRES CONFIGURABLES ==========

% Organisation par arbre (pour vérification)
trees = struct(...
    'Arbre1', {{'Arbre1_Nord', 'Arbre1_Est'}}, ...
    'Arbre2', {{'Arbre2_Nord', 'Arbre2_Est'}} ...
);

% ================================================

%% 1. Sélection du dossier contenant les fichiers MAT individuels

fprintf('\n=== Sélection du dossier contenant les fichiers MAT downsamplés ===\n');
input_folder = uigetdir('', 'Sélectionnez le dossier contenant les fichiers MAT individuels');

if input_folder == 0
    error('Opération annulée par l''utilisateur.');
end

fprintf('Dossier source : %s\n', input_folder);

%% 2. Sélection du dossier de sortie

fprintf('\n=== Sélection du dossier de sortie ===\n');
output_folder = uigetdir(input_folder, 'Sélectionnez le dossier pour enregistrer les fichiers MAT regroupés');

if output_folder == 0
    error('Opération annulée par l''utilisateur.');
end

fprintf('Dossier de sortie : %s\n', output_folder);

%% 3. Scan des fichiers MAT

fprintf('\n=== Scan des fichiers MAT ===\n');
mat_files = dir(fullfile(input_folder, '*.mat'));

if isempty(mat_files)
    error('Aucun fichier MAT trouvé dans le dossier %s', input_folder);
end

fprintf('Nombre de fichiers MAT trouvés : %d\n', length(mat_files));

%% 4. Groupement des fichiers par date/heure

fprintf('\n=== Groupement des fichiers par date/heure ===\n');

% Structure pour stocker les groupes : clé = date_str, valeur = liste de fichiers
file_groups = containers.Map('KeyType', 'char', 'ValueType', 'any');

for i = 1:length(mat_files)
    file_name = mat_files(i).name;
    
    % Extraction de la date/heure du nom de fichier
    % Format attendu : yyyymmdd_HHMMSS_ChannelName.mat
    tokens = regexp(file_name, '(\d{8}_\d{6})_(.+)\.mat', 'tokens');
    
    if isempty(tokens)
        fprintf('  ATTENTION : Format non reconnu pour %s\n', file_name);
        continue;
    end
    
    date_str = tokens{1}{1};  % yyyymmdd_HHMMSS
    channel_name = tokens{1}{2};  % Nom du canal
    
    % Ajouter au groupe correspondant
    if ~isKey(file_groups, date_str)
        file_groups(date_str) = struct('files', {}, 'channels', {});
    end
    
    current_group = file_groups(date_str);
    current_group(end+1).file_path = fullfile(mat_files(i).folder, file_name);
    current_group(end).channel_name = channel_name;
    file_groups(date_str) = current_group;
end

date_groups = keys(file_groups);
fprintf('Nombre de groupes date/heure trouvés : %d\n', length(date_groups));

%% 5. Concaténation et sauvegarde

fprintf('\n=== Concaténation des fichiers par date/heure ===\n');

for i = 1:length(date_groups)
    date_str = date_groups{i};
    group_files = file_groups(date_str);
    
    fprintf('[%d/%d] Traitement du groupe : %s\n', i, length(date_groups), date_str);
    
    % Structure pour stocker les données de tous les canaux
    data_struct = struct();
    date_heure = [];
    fs = [];
    
    % Charger tous les fichiers du groupe
    for j = 1:length(group_files)
        file_path = group_files(j).file_path;
        channel_name = group_files(j).channel_name;
        
        try
            % Charger le fichier MAT
            loaded_data = load(file_path);
            
            % Vérifier que les variables attendues existent
            if ~isfield(loaded_data, 'signal')
                fprintf('  ATTENTION : Variable "signal" manquante dans %s\n', file_path);
                continue;
            end
            
            % Stocker le signal avec le nom du canal
            data_struct.(channel_name) = loaded_data.signal;
            
            % Récupérer date_heure et fs (identiques pour tous les canaux du groupe)
            if isempty(date_heure) && isfield(loaded_data, 'date_heure')
                date_heure = loaded_data.date_heure;
            end
            if isempty(fs) && isfield(loaded_data, 'fs')
                fs = loaded_data.fs;
            end
            
            fprintf('  -> Canal %s : %d points chargés\n', channel_name, length(loaded_data.signal));
            
        catch ME
            fprintf('  ERREUR lors du chargement de %s : %s\n', file_path, ME.message);
            continue;
        end
    end
    
    % Vérifier qu'on a bien des données à sauvegarder
    if isempty(fieldnames(data_struct))
        fprintf('  -> Aucune donnée valide pour ce groupe, ignoré\n');
        continue;
    end
    
    % Créer le fichier de sortie
    output_filename = sprintf('%s.mat', date_str);
    output_path = fullfile(output_folder, output_filename);
    
    % Sauvegarder toutes les données dans un seul fichier
    % Structure finale : chaque canal est une variable + date_heure + fs
    save_vars = fieldnames(data_struct);
    for k = 1:length(save_vars)
        eval(sprintf('%s = data_struct.%s;', save_vars{k}, save_vars{k}));
    end
    
    % Sauvegarder avec toutes les variables
    save(output_path, save_vars{:}, 'date_heure', 'fs');
    
    fprintf('  => Fichier sauvegardé : %s (%d canaux)\n', output_filename, length(save_vars));
    fprintf('     Canaux inclus : %s\n', strjoin(save_vars, ', '));
end

fprintf('\n=== Concaténation terminée ! ===\n');
fprintf('Fichiers regroupés sauvegardés dans : %s\n', output_folder);

%% 6. Résumé final

fprintf('\n=== RÉSUMÉ ===\n');
fprintf('Total de groupes date/heure traités : %d\n', length(date_groups));
fprintf('Dossier de sortie : %s\n', output_folder);
fprintf('\nChaque fichier MAT contient :\n');
fprintf('  - Les signaux de tous les canaux (Arbre1_Nord, Arbre1_Est, Arbre2_Nord, Arbre2_Est)\n');
fprintf('  - La date/heure (variable "date_heure")\n');
fprintf('  - La fréquence d''échantillonnage (variable "fs")\n');