% ============================================================
% export_fmax_results.m
% ============================================================
% Sauvegarde les vecteurs (dates - fmax) pour chaque arbre et direction
% dans un fichier .mat compact et directement réutilisable.
%
% ⚙️ Ce script doit être exécuté APRÈS le calcul de 'results_daily'
%     dans le script principal d’analyse TDMS.
%
% 📦 Le fichier généré contiendra une structure :
%     export_data.Arbre1.Nord.dates
%     export_data.Arbre1.Nord.fmax
%     export_data.Arbre1.Est.dates
%     export_data.Arbre1.Est.fmax
%     etc.
%
% 🧩 Pour recharger plus tard :
%     load('V:\MONITORING_ARBRES\RESULTATS_FMAX\fmax_results.mat');
%     % Puis accéder ainsi :
%     dates = export_data.Arbre1.Nord.dates;
%     fmax  = export_data.Arbre1.Nord.fmax;
%
% Auteur : Lino CONORD
% Date   : 2025-10-15
% ============================================================


%% ==== PARAMÈTRES ====
output_folder = "V:\MONITORING_ARBRES\RESULTATS_FMAX";  % dossier de sortie
output_filename = fullfile(output_folder, "fmax_results.mat");

% Vérifie la présence de la variable 'results_daily' dans l’espace de travail
if ~evalin('base', 'exist(''results_daily'', ''var'')')
    error("La variable 'results_daily' n'existe pas dans l'espace de travail MATLAB.");
end

% Récupération depuis l’espace de travail principal
results_daily = evalin('base', 'results_daily');
tree_names = fieldnames(results_daily);

%% ==== Construction de la structure exportable ====
fprintf('\n=== Construction des vecteurs à exporter ===\n');

export_data = struct();

for ti = 1:length(tree_names)
    tree_name = tree_names{ti};
    directions = fieldnames(results_daily.(tree_name));
    
    for di = 1:length(directions)
        dir_name = directions{di};
        entry = results_daily.(tree_name).(dir_name);
        
        if isempty(entry.dates)
            fprintf('  %s - %s : aucun point à exporter.\n', tree_name, dir_name);
            continue;
        end
        
        % Conversion en struct de vecteurs simples
        export_data.(tree_name).(dir_name).dates = entry.dates(:);
        export_data.(tree_name).(dir_name).fmax  = entry.fmax(:);
        
        fprintf('  %s - %s : %d points exportés.\n', ...
            tree_name, dir_name, length(entry.dates));
    end
end

%% ==== Sauvegarde du fichier .mat ====
if ~isfolder(output_folder)
    mkdir(output_folder);
end

save(output_filename, 'export_data');

fprintf('\n✅ Fichier sauvegardé : %s\n', output_filename);
fprintf('Contenu : export_data.(Arbre).(Direction).dates / fmax\n');
