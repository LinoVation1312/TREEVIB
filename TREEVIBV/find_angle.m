% Lecture et traitement des fichiers MAT temporels avec visualisation interactive et calcul de direction
clear all; close all; clc;

%% ========== PARAMÈTRES CONFIGURABLES ==========
mat_folder = "C:\Users\lconord.LYD\OneDrive - MCO\Documents\MATLAB\treevib\data de Mars à Novembre 2025 - fs=25Hz\data de Mars à Novembre 2025 - fs=25Hz"; 
% Note : Modifiez le chemin ci-dessus si nécessaire

% Paramètres du signal
fs_target = 25;             
default_fs = 1000;          
f_dc = 0.12;                
fc = 0.6;                   
hp_taps = 2001;             
lp_order = 2;               

% Noms des canaux à lire (Structure pour mapping)
trees = struct(...
    'Arbre1', {{'Arbre1_Nord', 'Arbre1_Est'}}, ...
    'Arbre2', {{'Arbre2_Nord', 'Arbre2_Est'}} ...
);
colors_dir = [0 0.447 0.741; 0.85 0.33 0.1]; % Bleu (Nord), Rouge (Est)

%% ================================================
%% Vérification du dossier
if ~isfolder(mat_folder)
    warning("Le dossier %s est introuvable. Création de données simulées.", mat_folder);
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
        mat_files = all_files;
    end
else
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
        if isnat(dt), dt = datetime('now'); end 
        
        fs0 = default_fs;
        if isfield(S,'fs'), fs0 = double(S.fs); elseif isfield(S,'FS'), fs0 = double(S.FS); end
        
        % Extraction et Filtrage
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
        % Génération de données simulées (Ellipse qui tourne)
        file_name_only = "Simulation";
        dt = datetime('now');
        fs_target = 25;
        t = (0:1/fs_target:300)'; 
        
        % Arbre 1 : Ellipse inclinée à 30 degrés
        angle_tilt = deg2rad(30); 
        traj_x = 0.5*sin(2*pi*0.25*t);
        traj_y = 0.2*cos(2*pi*0.25*t);
        % Rotation
        A1_E_bp = traj_x * cos(angle_tilt) - traj_y * sin(angle_tilt) + 0.05*randn(size(t));
        A1_N_bp = traj_x * sin(angle_tilt) + traj_y * cos(angle_tilt) + 0.05*randn(size(t));
        
        % Arbre 2 : Ellipse inclinée à -45 degrés
        A2_N_bp = 0.4*sin(2*pi*0.3*t);
        A2_E_bp = 0.4*sin(2*pi*0.3*t + pi/2 + deg2rad(20)); % Phase shift complexe
    end
    
    %% === Lancement de l'interface interactive ===
    create_interactive_viewer(t, A1_N_bp, A1_E_bp, ...
        sprintf('Arbre 1 - %s (%s)', file_name_only, datestr(dt)), colors_dir);
        
    create_interactive_viewer(t, A2_N_bp, A2_E_bp, ...
        sprintf('Arbre 2 - %s (%s)', file_name_only, datestr(dt)), colors_dir);
        
end
fprintf('\n=== Traitement terminé ===\n');


%% ========== FONCTIONS LOCALES ==========

function [t, data_bp] = process_channel(raw_data, fs0, fs_target, f_dc, fc, hp_taps, lp_order)
    raw_data = raw_data(:);
    N_orig = length(raw_data);
    N_targ = round(N_orig * fs_target / fs0);
    
    data_resampled = resample(raw_data, N_targ, N_orig);
    
    nyq = fs_target/2;
    if f_dc > 0
        Wn_hp = f_dc/nyq;
        b_hp = fir1(hp_taps-1, Wn_hp, 'high', hamming(hp_taps), 'scale');
        y_hp = filtfilt(b_hp, 1, data_resampled);
    else
        y_hp = data_resampled;
    end
    
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


%% ========== VISUALISATION INTERACTIVE AVEC PCA ==========
%% ========== VISUALISATION INTERACTIVE AVEC VARIANCE (PCA) ==========
function create_interactive_viewer(t, dataN, dataE, title_str, colors)
    fig = figure('Name', title_str, 'NumberTitle', 'off', 'Position', [100, 100, 1400, 700]);
    
    % --- 1. Subplot Temporel (Gauche) ---
    ax_time = subplot(1,2,1, 'Parent', fig);
    hold(ax_time, 'on'); grid(ax_time, 'on');
    
    pN = plot(ax_time, t, dataN, 'Color', colors(1,:), 'DisplayName', 'Nord');
    pE = plot(ax_time, t, dataE, 'Color', colors(2,:), 'DisplayName', 'Est');
    
    ylabel(ax_time, 'Amplitude [m/s^2]'); xlabel(ax_time, 'Temps [s]');
    title(ax_time, {title_str; 'Série Temporelle'}, 'Interpreter', 'none', 'FontSize', 10);
    legend(ax_time, [pN, pE], 'Nord', 'Est', 'Location', 'northeast');
    axis(ax_time, 'tight');
    
    xl_min = xline(ax_time, t(1), '--k', 't_{min}', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    xl_max = xline(ax_time, t(end), '--k', 't_{max}', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

    % --- 2. Subplot Orbital (Droite) ---
    ax_orb = subplot(1,2,2, 'Parent', fig);
    hold(ax_orb, 'on'); grid(ax_orb, 'on'); axis(ax_orb, 'equal');
    
    plot(ax_orb, dataE, dataN, '.', 'Color', [0.9 0.9 0.9], 'MarkerSize', 4, 'DisplayName', 'Total');
    p_orb = plot(ax_orb, dataE, dataN, 'b.-', 'LineWidth', 0.5, 'DisplayName', 'Sélection');
    
    % --- Éléments graphiques PCA (Variance) ---
    % Axe Principal (Rouge)
    p_major = plot(ax_orb, [0 0], [0 0], 'r-', 'LineWidth', 2, 'DisplayName', 'Axe Majeur (Dir.)');
    % Axe Secondaire (Vert - Variance transverse)
    p_minor = plot(ax_orb, [0 0], [0 0], 'g-', 'LineWidth', 2, 'DisplayName', 'Axe Mineur');
    
    % Texte Info
    t_info = text(ax_orb, 0, 0, '', 'Color', 'k', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5);
    
    xlabel(ax_orb, 'Est [m/s^2]'); ylabel(ax_orb, 'Nord [m/s^2]');
    title(ax_orb, 'Trajectoire & Variances (PCA)');
    legend(ax_orb, [p_major, p_minor], 'Location', 'southoutside', 'Orientation', 'horizontal');
    
    % --- Contrôles UI ---
    bg_color = [0.95 0.95 0.95];
    panel_height = 0.15;
    ui_panel = uipanel(fig, 'Position', [0 0 1 panel_height], 'BackgroundColor', bg_color);
    
    set(ax_time, 'Position', [0.05, 0.22, 0.42, 0.70]); 
    set(ax_orb,  'Position', [0.53, 0.22, 0.42, 0.70]);
    
    uicontrol(ui_panel, 'Style', 'text', 'Position', [20 50 50 20], 'String', 'T min:', 'BackgroundColor', bg_color);
    uicontrol(ui_panel, 'Style', 'text', 'Position', [20 20 50 20], 'String', 'T max:', 'BackgroundColor', bg_color);
    
    sl_min = uicontrol(ui_panel, 'Style', 'slider', 'Position', [80 50 600 20], ...
        'Min', t(1), 'Max', t(end), 'Value', t(1), 'Callback', @update_plot);
    sl_max = uicontrol(ui_panel, 'Style', 'slider', 'Position', [80 20 600 20], ...
        'Min', t(1), 'Max', t(end), 'Value', t(end), 'Callback', @update_plot);
        
    ed_min = uicontrol(ui_panel, 'Style', 'edit', 'Position', [700 50 80 20], ...
        'String', num2str(t(1)), 'Callback', @update_edit_min);
    ed_max = uicontrol(ui_panel, 'Style', 'edit', 'Position', [700 20 80 20], ...
        'String', num2str(t(end)), 'Callback', @update_edit_max);

    % --- Callback ---
    function update_plot(~, ~)
        val_min = get(sl_min, 'Value');
        val_max = get(sl_max, 'Value');
        if val_min >= val_max, val_min = val_max - 0.01; set(sl_min, 'Value', val_min); end
        
        set(ed_min, 'String', sprintf('%.2f', val_min));
        set(ed_max, 'String', sprintf('%.2f', val_max));
        xl_min.Value = val_min; xl_max.Value = val_max;
        
        idx = t >= val_min & t <= val_max;
        
        if sum(idx) > 2 % Besoin d'au moins 2 points pour covariance
            subE = dataE(idx); subN = dataN(idx);
            set(p_orb, 'XData', subE, 'YData', subN);
            
            % --- Calcul PCA ---
            mean_E = mean(subE); mean_N = mean(subN);
            C = cov(subE - mean_E, subN - mean_N);
            [V, D] = eig(C);
            
            % Tri des valeurs propres (du plus grand au plus petit)
            [eig_vals, sort_idx] = sort(diag(D), 'descend');
            eig_vecs = V(:, sort_idx);
            
            % Écart-type (racine de la variance) pour échelle visuelle
            % On multiplie par 2 pour couvrir ~95% des données (si gaussien)
            sigma_major = 2 * sqrt(eig_vals(1)); 
            sigma_minor = 2 * sqrt(eig_vals(2));
            
            % Vecteurs directeurs
            v_major = eig_vecs(:, 1);
            v_minor = eig_vecs(:, 2);
            
            % Calcul Angle Principal
            angle_deg = atan2d(v_major(2), v_major(1));
            % Normalisation -90 à 90
            if angle_deg > 90, angle_deg = angle_deg - 180; end
            if angle_deg < -90, angle_deg = angle_deg + 180; end
            
            % Ratio d'aspect (0 = ligne, 1 = cercle)
            ratio_aspect = sqrt(eig_vals(2) / eig_vals(1));
            
            % --- Mise à jour Graphique ---
            % Ligne Majeure (Rouge)
            set(p_major, 'XData', [mean_E - v_major(1)*sigma_major, mean_E + v_major(1)*sigma_major], ...
                         'YData', [mean_N - v_major(2)*sigma_major, mean_N + v_major(2)*sigma_major]);
            
            % Ligne Mineure (Verte)
            set(p_minor, 'XData', [mean_E - v_minor(1)*sigma_minor, mean_E + v_minor(1)*sigma_minor], ...
                         'YData', [mean_N - v_minor(2)*sigma_minor, mean_N + v_minor(2)*sigma_minor]);
            
            % Info Texte
            max_range = max(max(subE)-min(subE), max(subN)-min(subN));
            scale_pos = max_range * 0.55;
            info_str = {sprintf('Angle: %.1f°', angle_deg), ...
                        sprintf('Ratio (Min/Maj): %.2f', ratio_aspect)};
            
            set(t_info, 'Position', [mean_E, mean_N + scale_pos], 'String', info_str);
            
            % Zoom Auto
            padding = max_range * 0.2;
            xlim(ax_orb, [mean_E - max_range/2 - padding, mean_E + max_range/2 + padding]);
            ylim(ax_orb, [mean_N - max_range/2 - padding, mean_N + max_range/2 + padding]);
        end
    end

    function update_edit_min(~, ~), val = str2double(get(ed_min,'String')); set(sl_min,'Value',min(max(val,t(1)),t(end))); update_plot(); end
    function update_edit_max(~, ~), val = str2double(get(ed_max,'String')); set(sl_max,'Value',min(max(val,t(1)),t(end))); update_plot(); end

    update_plot();
end