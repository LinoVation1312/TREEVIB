% Lecture et traitement des fichiers MAT temporels avec visualisation interactive
clear all; close all; clc;

%% ========== PARAMÈTRES CONFIGURABLES ==========
mat_folder = "H:\Home\Documents\aaaTREEVIB\downsampling_filtered"; 
% Note : Modifiez le chemin ci-dessus si nécessaire

% Paramètres du signal
fs_target = 25;             
default_fs = 1000;          
f_dc = 0.12;                
fc = 0.6;                   
hp_taps = 2001;             
lp_order = 2;               

% Paramètres de pwelch (pour l'analyse spectrale en arrière-plan si besoin)
pwelch_window = 300 * fs_target;
pwelch_overlap = 0;
pwelch_nfft = 2^nextpow2(pwelch_window);

% Noms des canaux à lire
trees = struct(...
    'Arbre1', {{'Arbre1_Nord', 'Arbre1_Est'}}, ...
    'Arbre2', {{'Arbre2_Nord', 'Arbre2_Est'}} ...
);
channel_names = unique([trees.Arbre1, trees.Arbre2]); 
colors_dir = [0 0.447 0.741; 0.85 0.33 0.1]; % Bleu (Nord), Rouge (Est)

%% ================================================
%% Vérification du dossier
if ~isfolder(mat_folder)
    % Fallback pour l'exemple si le dossier n'existe pas chez l'utilisateur
    warning("Le dossier %s est introuvable. Création de données simulées pour démo.", mat_folder);
    USE_SIMULATION = true;
else
    USE_SIMULATION = false;
end

%% Scan et Sélection des fichiers MAT
if ~USE_SIMULATION
    fprintf('\n=== Scan des fichiers MAT dans le dossier ===\n');
    all_files = dir(fullfile(mat_folder, '*.mat'));
    
    if isempty(all_files)
        error('Aucun fichier MAT trouvé dans le dossier %s', mat_folder);
    end
    
    [~, sort_idx] = sort({all_files.name});
    all_files = all_files(sort_idx);
    file_names_list = {all_files.name}';
    
    try
        [selection_indices, ok] = listdlg('PromptString','Sélectionner les fichiers MAT à analyser:', ...
                                          'SelectionMode','multiple', ...
                                          'ListString', file_names_list, ...
                                          'Name','Sélection des Données');
        if ~ok || isempty(selection_indices), error('Aucune sélection.'); end
        mat_files = all_files(selection_indices);
    catch
        mat_files = all_files; % Fallback
    end
else
    % Création d'une structure factice pour la simulation sans fichiers
    mat_files = struct('name', 'Simulation_Data.mat', 'folder', pwd);
end

%% ========== BOUCLE PRINCIPALE DE TRAITEMENT ==========
for k = 1:numel(mat_files)
    if ~USE_SIMULATION
        fname = fullfile(mat_files(k).folder, mat_files(k).name);
        file_name_only = mat_files(k).name;
        fprintf('\n[%d/%d] Traitement de %s\n', k, numel(mat_files), file_name_only);
        
        try
            S = load(fname);
        catch ME
            warning('Erreur chargement %s : %s', file_name_only, ME.message);
            continue;
        end
        
        % Date et Fs
        dt = parse_filename_date(file_name_only, S);
        if isnat(dt), dt = datetime('now'); end % Fallback date
        
        fs0 = default_fs;
        if isfield(S,'fs'), fs0 = double(S.fs); elseif isfield(S,'FS'), fs0 = double(S.FS); end
        
        % Extraction et Filtrage des 4 canaux
        % (On suppose ici que les 4 canaux existent pour simplifier le code interactif)
        try
            % Arbre 1
            [t, A1_N_bp] = process_channel(S.Arbre1_Nord, fs0, fs_target, f_dc, fc, hp_taps, lp_order);
            [~, A1_E_bp] = process_channel(S.Arbre1_Est, fs0, fs_target, f_dc, fc, hp_taps, lp_order);
            
            % Arbre 2
            [~, A2_N_bp] = process_channel(S.Arbre2_Nord, fs0, fs_target, f_dc, fc, hp_taps, lp_order);
            [~, A2_E_bp] = process_channel(S.Arbre2_Est, fs0, fs_target, f_dc, fc, hp_taps, lp_order);
        catch
            warning('Canaux manquants dans %s. Passage au suivant.', file_name_only);
            continue;
        end
        
    else
        % Génération de données simulées
        file_name_only = "Simulation";
        dt = datetime('now');
        fs_target = 25;
        t = (0:1/fs_target:300)'; % 5 minutes
        
        % Arbre 1 (Mouvement elliptique simple + bruit)
        A1_N_bp = 0.5*sin(2*pi*0.25*t) + 0.1*randn(size(t));
        A1_E_bp = 0.3*cos(2*pi*0.25*t + pi/4) + 0.1*randn(size(t));
        
        % Arbre 2 (Mouvement plus chaotique)
        A2_N_bp = 0.4*sin(2*pi*0.3*t);
        A2_E_bp = 0.4*sin(2*pi*0.35*t);
    end
    
    %% === Lancement de l'interface interactive ===
    % On crée une figure interactive pour Arbre 1
    create_interactive_viewer(t, A1_N_bp, A1_E_bp, ...
        sprintf('Arbre 1 - %s (%s)', file_name_only, datestr(dt)), colors_dir);
        
    % On crée une figure interactive pour Arbre 2
    create_interactive_viewer(t, A2_N_bp, A2_E_bp, ...
        sprintf('Arbre 2 - %s (%s)', file_name_only, datestr(dt)), colors_dir);
        
end
fprintf('\n=== Traitement terminé ===\n');


%% ========== FONCTIONS LOCALES ==========

function [t, data_bp] = process_channel(raw_data, fs0, fs_target, f_dc, fc, hp_taps, lp_order)
    % Resampling et Filtrage
    raw_data = raw_data(:);
    N_orig = length(raw_data);
    N_targ = round(N_orig * fs_target / fs0);
    
    data_resampled = resample(raw_data, N_targ, N_orig);
    
    % Filtrage (Copie de la logique précédente)
    nyq = fs_target/2;
    
    % HP
    if f_dc > 0
        Wn_hp = f_dc/nyq;
        b_hp = fir1(hp_taps-1, Wn_hp, 'high', hamming(hp_taps), 'scale');
        y_hp = filtfilt(b_hp, 1, data_resampled);
    else
        y_hp = data_resampled;
    end
    
    % LP
    Wn_lp = fc/nyq;
    [b_lp, a_lp] = butter(lp_order, Wn_lp, 'low');
    data_bp = filtfilt(b_lp, a_lp, y_hp);
    
    t = (0:length(data_bp)-1)' / fs_target;
end

function dt = parse_filename_date(fname, S)
    dt = NaT;
    if isfield(S,'date_heure'), try dt = datetime(S.date_heure); catch, end; end
    if isnat(dt)
        tok = regexp(fname, '(\d{4})(\d{2})(\d{2})[_-]?(\d{2})(\d{2})(\d{2})', 'tokens');
        if ~isempty(tok), v=tok{1}; dt = datetime(str2double(v{1}),str2double(v{2}),str2double(v{3}),str2double(v{4}),str2double(v{5}),str2double(v{6})); end
    end
end


%% ========== VISUALISATION INTERACTIVE ==========
function create_interactive_viewer(t, dataN, dataE, title_str, colors)
    % Création de la figure
    % On élargit un peu la figure par défaut pour le layout côte à côte
    fig = figure('Name', title_str, 'NumberTitle', 'off', 'Position', [100, 100, 1400, 700]);
    
    % --- Layout Graphique : Côte à Côte ---
    
    % 1. Subplot Temporel (Gauche) -> 1,2,1
    ax_time = subplot(1,2,1, 'Parent', fig);
    hold(ax_time, 'on'); grid(ax_time, 'on');
    
    % On garde les handles pN et pE pour la légende
    pN = plot(ax_time, t, dataN, 'Color', colors(1,:), 'DisplayName', 'Nord');
    pE = plot(ax_time, t, dataE, 'Color', colors(2,:), 'DisplayName', 'Est');
    
    ylabel(ax_time, 'Amplitude [m/s^2]');
    xlabel(ax_time, 'Temps [s]');
    
    % Titre avec le nom de l'arbre
    title(ax_time, {title_str; 'Signaux Temporels Filtrés'}, 'Interpreter', 'none', 'FontSize', 10);
    
    % CORRECTION ICI : On force la légende uniquement sur pN et pE
    legend(ax_time, [pN, pE], 'Nord', 'Est', 'Location', 'northeast');
    
    axis(ax_time, 'tight');
    
    % Lignes verticales pour curseurs
    % On ajoute 'HandleVisibility','off' pour être sûr qu'elles n'apparaissent jamais dans la légende
    t_min_init = t(1);
    t_max_init = t(end);
    xl_min = xline(ax_time, t_min_init, '--k', 't_{min}', ...
        'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    xl_max = xline(ax_time, t_max_init, '--k', 't_{max}', ...
        'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

    % 2. Subplot Orbital (Droite) -> 1,2,2
    ax_orb = subplot(1,2,2, 'Parent', fig);
    
    % On trace tout en gris clair pour référence
    plot(ax_orb, dataE, dataN, '.', 'Color', [0.8 0.8 0.8], 'MarkerSize', 4, 'DisplayName', 'Total');
    hold(ax_orb, 'on'); grid(ax_orb, 'on');
    p_orb = plot(ax_orb, dataE, dataN, 'b.-', 'LineWidth', 0.5, 'DisplayName', 'Sélection');
    
    xlabel(ax_orb, 'Est [m/s^2]');
    ylabel(ax_orb, 'Nord [m/s^2]');
    
    % Titre avec le nom de l'arbre
    title(ax_orb, {title_str; 'Plan Horizontal (Vue de dessus)'}, 'Interpreter', 'none', 'FontSize', 10);
    
    axis(ax_orb, 'equal'); % MÊME ÉCHELLE X et Y
    
    % --- Contrôles UI (Sliders) ---
    bg_color = [0.95 0.95 0.95];
    panel_height = 0.15;
    ui_panel = uipanel(fig, 'Position', [0 0 1 panel_height], 'BackgroundColor', bg_color);
    
    % Ajustement manuel des positions des axes
    set(ax_time, 'Position', [0.05, 0.22, 0.42, 0.70]); 
    set(ax_orb,  'Position', [0.53, 0.22, 0.42, 0.70]);
    
    % Labels
    uicontrol(ui_panel, 'Style', 'text', 'Position', [20 50 50 20], 'String', 'T min:', 'BackgroundColor', bg_color);
    uicontrol(ui_panel, 'Style', 'text', 'Position', [20 20 50 20], 'String', 'T max:', 'BackgroundColor', bg_color);
    
    % Sliders
    sl_min = uicontrol(ui_panel, 'Style', 'slider', 'Position', [80 50 600 20], ...
        'Min', t(1), 'Max', t(end), 'Value', t(1), 'Callback', @update_plot);
    sl_max = uicontrol(ui_panel, 'Style', 'slider', 'Position', [80 20 600 20], ...
        'Min', t(1), 'Max', t(end), 'Value', t(end), 'Callback', @update_plot);
        
    % Edits
    ed_min = uicontrol(ui_panel, 'Style', 'edit', 'Position', [700 50 80 20], ...
        'String', num2str(t(1)), 'Callback', @update_edit_min);
    ed_max = uicontrol(ui_panel, 'Style', 'edit', 'Position', [700 20 80 20], ...
        'String', num2str(t(end)), 'Callback', @update_edit_max);

    % --- Fonctions de Callback (Nested) ---
    function update_plot(~, ~)
        val_min = get(sl_min, 'Value');
        val_max = get(sl_max, 'Value');
        
        if val_min >= val_max
            val_min = val_max - (t(2)-t(1)); 
            set(sl_min, 'Value', val_min);
        end
        
        set(ed_min, 'String', sprintf('%.2f', val_min));
        set(ed_max, 'String', sprintf('%.2f', val_max));
        
        xl_min.Value = val_min;
        xl_max.Value = val_max;
        
        idx = t >= val_min & t <= val_max;
        
        if any(idx)
            set(p_orb, 'XData', dataE(idx), 'YData', dataN(idx));
            
            x_sel = dataE(idx); y_sel = dataN(idx);
            max_range = max(max(x_sel)-min(x_sel), max(y_sel)-min(y_sel));
            mid_x = (max(x_sel)+min(x_sel))/2;
            mid_y = (max(y_sel)+min(y_sel))/2;
            
            padding = max_range * 0.1; 
            if padding == 0, padding = 0.1; end
            half_span = (max_range/2) + padding;
            
            xlim(ax_orb, [mid_x - half_span, mid_x + half_span]);
            ylim(ax_orb, [mid_y - half_span, mid_y + half_span]);
        end
    end

    function update_edit_min(~, ~)
        val = str2double(get(ed_min, 'String'));
        if val < t(1), val = t(1); end
        if val > get(sl_max, 'Value'), val = get(sl_max, 'Value') - 0.1; end
        set(sl_min, 'Value', val);
        update_plot();
    end

    function update_edit_max(~, ~)
        val = str2double(get(ed_max, 'String'));
        if val > t(end), val = t(end); end
        if val < get(sl_min, 'Value'), val = get(sl_min, 'Value') + 0.1; end
        set(sl_max, 'Value', val);
        update_plot();
    end

    % Initial call
    update_plot();
end