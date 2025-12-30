% Analyse de tendance : Angle de vibration vs Temps pour tout le dossier
clear all; close all; clc;

%% ========== 1. PARAMÈTRES ==========
mat_folder = "C:\Users\lconord.LYD\OneDrive - MCO\Documents\MATLAB\treevib\data de Mars à Novembre 2025 - fs=25Hz\data de Mars à Novembre 2025 - fs=25Hz"; 

% Paramètres de traitement
fs_target = 25;             
default_fs = 1000;          
f_dc = 0.12;                
fc = 0.6;                   
hp_taps = 2001;             
lp_order = 2;               

%% ========== 2. PRÉPARATION ==========
if ~isfolder(mat_folder)
    % Mode simulation si dossier absent
    warning('Dossier introuvable. Mode simulation.');
    files = struct('name', 'Simu.mat', 'folder', pwd);
    IS_SIMU = true;
else
    files = dir(fullfile(mat_folder, '*.mat'));
    IS_SIMU = false;
end

if isempty(files) && ~IS_SIMU
    error('Aucun fichier .mat trouvé.');
end

results_A1 = []; 
results_A2 = [];

fprintf('Début du traitement...\n');
wb = waitbar(0, 'Traitement en cours...');

%% ========== 3. BOUCLE DE TRAITEMENT ==========
for k = 1:length(files)
    waitbar(k/length(files), wb, sprintf('Fichier %d/%d', k, length(files)));
    
    if ~IS_SIMU
        fname = files(k).name;
        full_path = fullfile(files(k).folder, fname);
        try
            S = load(full_path);
        catch
            continue; 
        end
        
        dt = parse_filename_date(fname, S);
        if isnat(dt), continue; end 
        t_num = datenum(dt);
        
        fs0 = default_fs;
        if isfield(S,'fs'), fs0 = double(S.fs); elseif isfield(S,'FS'), fs0 = double(S.FS); end
        
        % Arbre 1
        if isfield(S, 'Arbre1_Nord') && isfield(S, 'Arbre1_Est')
            [res_angle, res_ratio] = compute_pca_angle(S.Arbre1_Nord, S.Arbre1_Est, fs0, fs_target, f_dc, fc, hp_taps, lp_order);
            % On ne garde que si l'angle est un nombre valide (pas de NaN)
            if ~isnan(res_angle)
                results_A1 = [results_A1; t_num, res_angle, res_ratio];
            end
        end
        
        % Arbre 2
        if isfield(S, 'Arbre2_Nord') && isfield(S, 'Arbre2_Est')
            [res_angle, res_ratio] = compute_pca_angle(S.Arbre2_Nord, S.Arbre2_Est, fs0, fs_target, f_dc, fc, hp_taps, lp_order);
            if ~isnan(res_angle)
                results_A2 = [results_A2; t_num, res_angle, res_ratio];
            end
        end
        
    else
        % Simulation de données pour l'exemple
        t_num = datenum(datetime('now') + hours(k));
        results_A1 = [results_A1; t_num, 45 + 5*randn, 0.6 + 0.1*randn];
        results_A2 = [results_A2; t_num, -30 + 5*randn, 0.4 + 0.1*randn];
    end
end
close(wb);

%% ========== 4. VISUALISATION ==========
fprintf('Génération des graphiques...\n');

% On trace seulement s'il y a des données
if ~isempty(results_A1)
    plot_trend(results_A1, 'Arbre 1', [0 0.447 0.741]);
else
    disp('Pas de données valides pour Arbre 1');
end

if ~isempty(results_A2)
    plot_trend(results_A2, 'Arbre 2', [0.85 0.33 0.1]);
else
    disp('Pas de données valides pour Arbre 2');
end


%% ========== FONCTIONS LOCALES ==========

function [angle_deg, ratio] = compute_pca_angle(raw_N, raw_E, fs0, fs_target, f_dc, fc, hp_taps, lp_order)
    try
        [~, N_bp] = process_channel(raw_N, fs0, fs_target, f_dc, fc, hp_taps, lp_order);
        [~, E_bp] = process_channel(raw_E, fs0, fs_target, f_dc, fc, hp_taps, lp_order);
        
        if isempty(N_bp) || isempty(E_bp)
            angle_deg = NaN; ratio = NaN; return; 
        end
        
        % PCA
        data = [E_bp(:) - mean(E_bp), N_bp(:) - mean(N_bp)];
        C = cov(data);
        [V, D] = eig(C);
        [vals, idx] = sort(diag(D), 'descend');
        vecs = V(:, idx);
        
        principal_vec = vecs(:, 1);
        angle_deg = atan2d(principal_vec(2), principal_vec(1));
        
        % Normalisation -90 à 90
        if angle_deg > 90, angle_deg = angle_deg - 180; end
        if angle_deg < -90, angle_deg = angle_deg + 180; end
        
        ratio = sqrt(vals(2) / vals(1));
        
    catch
        angle_deg = NaN; ratio = NaN;
    end
end

function [t, data_bp] = process_channel(raw_data, fs0, fs_target, f_dc, fc, hp_taps, lp_order)
    raw_data = raw_data(:);
    N_orig = length(raw_data);
    N_targ = round(N_orig * fs_target / fs0);
    if N_targ < 10, t=[]; data_bp=[]; return; end % Protection fichier vide
    
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

function plot_trend(results, tree_name, main_color)
    % Tri par date
    [~, idx] = sort(results(:,1));
    results = results(idx, :);
    
    dates_num = results(:,1);
    angles = results(:,2);
    ratios = results(:,3);
    
    figure('Name', ['Tendance - ' tree_name], 'NumberTitle', 'off', 'Position', [100 100 1000 600]);
    
    % --- Subplot 1 : Angle avec couleur selon Ratio ---
    subplot(2,1,1);
    
    % Scatter : Couleur bornée
    scatter(datetime(dates_num, 'ConvertFrom', 'datenum'), angles, 30, ratios, 'filled');
    
    % CONFIGURATION COULEUR DEMANDÉE
    colormap(flipud(parula)); 
    clim([0.5 1]); % Borne stricte entre 0.5 et 1
    cb = colorbar; 
    ylabel(cb, 'Ratio (0.5=Plat/Fiable, 1=Rond)');
    
    hold on;
    
    % --- Ligne de tendance avec gestion des TROUS ---
    % Si l'écart entre deux points est > 2x la médiane, on coupe la ligne (NaN)
    if length(dates_num) > 1
        dt_diff = diff(dates_num);
        median_diff = median(dt_diff);
        
        % On insère des NaNs là où il y a des trous pour couper le plot
        dates_clean = dates_num;
        angles_clean = movmean(angles, 5); % Lissage léger
        
        gap_idx = find(dt_diff > 1.5 * median_diff);
        
        % Insertion astucieuse de NaNs pour le plot
        plot_t = [];
        plot_a = [];
        last_idx = 1;
        
        for i = 1:length(gap_idx)
            current_gap = gap_idx(i);
            % Ajouter le segment valide
            plot_t = [plot_t; dates_clean(last_idx:current_gap); NaN];
            plot_a = [plot_a; angles_clean(last_idx:current_gap); NaN];
            last_idx = current_gap + 1;
        end
        % Dernier segment
        plot_t = [plot_t; dates_clean(last_idx:end)];
        plot_a = [plot_a; angles_clean(last_idx:end)];
        
        plot(datetime(plot_t, 'ConvertFrom', 'datenum'), plot_a, 'k-', 'LineWidth', 1.5);
    else
        plot(datetime(dates_num, 'ConvertFrom', 'datenum'), angles, 'k-');
    end
    
    ylabel('Angle [Deg]');
    title(['Direction - ' tree_name]);
    grid on; ylim([-95 95]); yline(0, 'k--');
    
    % --- Subplot 2 : Ratio ---
    subplot(2,1,2);
    plot(datetime(dates_num, 'ConvertFrom', 'datenum'), ratios, '.-', 'Color', main_color);
    ylabel('Ratio');
    xlabel('Temps');
    grid on; ylim([0 1]);
end