% Histogramme de disponibilité des données TDMS par jour
clear all; close all; clc;

%% ========== PARAMÈTRES CONFIGURABLES ==========
tdms_folder = "V:\MONITORING_ARBRES\PEUPLIERS a partir du 03 mars (a mettre a jour chaque mois)";

%% Vérification du dossier
if ~isfolder(tdms_folder)
    error("Le dossier %s est introuvable.", tdms_folder);
end

%% Scan des fichiers TDMS
fprintf('\n=== Scan des fichiers TDMS ===\n');
tdms_files = dir(fullfile(tdms_folder, '*.tdms'));

if isempty(tdms_files)
    error('Aucun fichier TDMS trouvé dans le dossier %s', tdms_folder);
end

fprintf('Nombre de fichiers TDMS trouvés : %d\n', length(tdms_files));

%% Extraction des dates depuis les noms de groupes
all_dates = [];
failed_count = 0;

fprintf('\n=== Extraction des dates depuis les groupes TDMS ===\n');

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
        
        % Récupérer tous les groupes uniques du fichier
        file_groups = unique(chanlist.ChannelGroupName);
        
        % Parser chaque nom de groupe pour extraire la date
        for g = 1:length(file_groups)
            group_name = file_groups{g};
            
            % Format attendu: YYYY_MM_DD_HH_MM_SS
            tokens = regexp(group_name, '(\d{4})_(\d{2})_(\d{2})_(\d{2})_(\d{2})_(\d{2})', 'tokens');
            
            if ~isempty(tokens)
                % Extraire année, mois, jour, heure, minute, seconde
                t = tokens{1};
                year = str2double(t{1});
                month = str2double(t{2});
                day = str2double(t{3});
                hour = str2double(t{4});
                minute = str2double(t{5});
                second = str2double(t{6});
                
                % Créer un datetime
                dt = datetime(year, month, day, hour, minute, second);
                all_dates = [all_dates; dt];
            else
                failed_count = failed_count + 1;
            end
        end
        
    catch ME
        fprintf('  Erreur fichier %s : %s\n', tdms_files(i).name, ME.message);
    end
end

if isempty(all_dates)
    error('Aucune date valide extraite des groupes TDMS');
end

fprintf('\nDates extraites : %d\n', length(all_dates));
fprintf('Groupes non parsés : %d\n', failed_count);
fprintf('Plage temporelle : %s à %s\n', datestr(min(all_dates)), datestr(max(all_dates)));

%% Agrégation par jour
% Convertir en dates (sans heure)
dates_only = dateshift(all_dates, 'start', 'day');
unique_days = unique(dates_only);

fprintf('\nNombre de jours distincts avec données : %d\n', length(unique_days));

%% Créer une série temporelle complète (tous les jours)
date_min = min(unique_days);
date_max = max(unique_days);
all_days = (date_min:days(1):date_max)';

% Vecteur de disponibilité (1 si données, 0 sinon)
availability = ismember(all_days, unique_days);

fprintf('Période totale : %d jours\n', length(all_days));
fprintf('Jours avec données : %d (%.1f%%)\n', sum(availability), 100*sum(availability)/length(all_days));
fprintf('Jours sans données : %d (%.1f%%)\n', sum(~availability), 100*sum(~availability)/length(all_days));

%% Visualisation
figure('Position', [100 100 1400 600]);

% Subplot 1: Histogramme de disponibilité
subplot(2,1,1);
bar(all_days, availability, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none', 'BarWidth', 1);
grid on;
xlabel('Date', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Disponibilité', 'FontSize', 12, 'FontWeight', 'bold');
title('Disponibilité des données TDMS par jour (1 = données présentes, 0 = absent)', ...
      'FontSize', 14, 'FontWeight', 'bold');
ylim([-0.1 1.1]);
yticks([0 1]);
yticklabels({'Absent', 'Présent'});

% Ajuster l'axe des dates
xlim([date_min-days(1), date_max+days(1)]);

% Subplot 2: Nombre de mesures par jour
measurements_per_day = histcounts(dates_only, [unique_days; unique_days(end)+days(1)]);
subplot(2,1,2);
bar(unique_days, measurements_per_day, 'FaceColor', [0.8 0.4 0.2], 'EdgeColor', 'none');
grid on;
xlabel('Date', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Nombre de mesures', 'FontSize', 12, 'FontWeight', 'bold');
title('Nombre de groupes TDMS (mesures) par jour', 'FontSize', 14, 'FontWeight', 'bold');
xlim([date_min-days(1), date_max+days(1)]);

%% Statistiques détaillées
fprintf('\n=== STATISTIQUES DÉTAILLÉES ===\n');
fprintf('Période d''analyse : %s à %s\n', datestr(date_min), datestr(date_max));
fprintf('Durée totale : %d jours\n', days(date_max - date_min) + 1);
fprintf('Taux de disponibilité : %.2f%%\n', 100*sum(availability)/length(all_days));

% Identifier les périodes sans données
availability_padded = [1; availability; 1];  % Pad avec 1 pour détecter correctement
gaps = diff(availability_padded);
gap_starts = find(gaps == -1);
gap_ends = find(gaps == 1) - 1;

if ~isempty(gap_starts) && ~isempty(gap_ends)
    fprintf('\n=== PÉRIODES SANS DONNÉES (lacunes) ===\n');
    for i = 1:min(length(gap_starts), length(gap_ends))
        if gap_starts(i) <= length(all_days) && gap_ends(i) <= length(all_days)
            gap_duration = gap_ends(i) - gap_starts(i) + 1;
            fprintf('  %s à %s (%d jours)\n', ...
                    datestr(all_days(gap_starts(i))), ...
                    datestr(all_days(gap_ends(i))), ...
                    gap_duration);
        end
    end
else
    fprintf('\nAucune lacune détectée : données continues!\n');
end

% Statistiques sur le nombre de mesures par jour
fprintf('\n=== STATISTIQUES DES MESURES PAR JOUR ===\n');
fprintf('Moyenne : %.1f mesures/jour\n', mean(measurements_per_day));
fprintf('Médiane : %.1f mesures/jour\n', median(measurements_per_day));
fprintf('Min : %d mesures/jour\n', min(measurements_per_day));
fprintf('Max : %d mesures/jour\n', max(measurements_per_day));
fprintf('Total : %d mesures\n', sum(measurements_per_day));