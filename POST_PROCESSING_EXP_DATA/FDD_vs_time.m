clearvars; close all; clc;

%% CONFIGURATION
MAT_FOLDER = "H:\Home\Documents\aaaTREEVIB\new_5Hz";

FS_RAW        = 5;    
DO_DOWNSAMPLE = false;  
FS_DS_VAL     = 5;    

if DO_DOWNSAMPLE
    FS_TARGET = FS_DS_VAL;
    fprintf("Mode: Downsampling activé (%d Hz -> %d Hz)\n", FS_RAW, FS_TARGET);
else
    FS_TARGET = FS_RAW;
    fprintf("Mode: Fréquence native conservée (%d Hz)\n", FS_TARGET);
end

TARGET_POINTS_RAW = 1500; % 300 sec × 5 Hz

% --- PARAMETRES DE L'IMAGE ---
DEFAULT_FMIN = 0.1;          
DEFAULT_FMAX = 0.75;         
DEFAULT_DELTAF = 0.1;        
DEFAULT_DYN = 10;            
DEFAULT_SVD_RATIO = 5;       
TARGET_NFFT_PADDING = 2048; 
% -----------------------------

BUTTER_ORDER = 2;
DISPLAY_FREQ_LIMIT = [0 1];  
Y_MARGIN_UP = 5;             
Y_DYNAMIC_RANGE = 30;   
GAP_THRESHOLD = hours(12);   
NUM_MODES_MAX = 5;         

%% INITIALISATION
TreeConfigs = struct();
TreeConfigs(1).Name = 'Arbre 1';
TreeConfigs(1).ChN  = 'Arbre1_Nord';
TreeConfigs(1).ChE  = 'Arbre1_Est';
TreeConfigs(2).Name = 'Arbre 2';
TreeConfigs(2).ChN  = 'Arbre2_Nord';
TreeConfigs(2).ChE  = 'Arbre2_Est';

files = dir(fullfile(MAT_FOLDER, "*.mat"));
if isempty(files), error("Aucun fichier .mat trouvé dans : %s", MAT_FOLDER); end

nFiles = numel(files);
idx_rand = randi(nFiles);

%% TRAITEMENT
for iTree = 1:length(TreeConfigs)
    cfg = TreeConfigs(iTree);
    fprintf("Traitement : %s ...\n", cfg.Name);
    
    all_dates = NaT(nFiles, 1);
    all_peaks = nan(nFiles, NUM_MODES_MAX); 
    snap = struct('valid', false);
    
    %ITERATION
    for k = 1:nFiles
        fname = files(k).name;
        if mod(k, 50) == 0, fprintf("  [%d/%d] \n", k, nFiles); end
        
        try
            S = load(fullfile(files(k).folder, fname));
        catch
            continue;
        end
        
        %DATE
        dt = NaT;
        if isfield(S,'date_heure'), try, dt = datetime(S.date_heure); catch, end; end
        if isnat(dt), dt = parse_filename_date(fname); end
        all_dates(k) = dt;
        
        %VERIFICATION
        if ~isfield(S, cfg.ChN) || ~isfield(S, cfg.ChE), continue; end
        
        [dN] = process_signal(S.(cfg.ChN), FS_RAW, FS_TARGET, DO_DOWNSAMPLE, BUTTER_ORDER, TARGET_POINTS_RAW);
        [dE] = process_signal(S.(cfg.ChE), FS_RAW, FS_TARGET, DO_DOWNSAMPLE, BUTTER_ORDER, TARGET_POINTS_RAW);
        
        nMin = min(length(dN), length(dE));
        dN = dN(1:nMin); dE = dE(1:nMin);
        
        % --- CONFIGURATION WELCH OPTIMALE POUR 1500 POINTS ---
        % Objectif : équilibre résolution fréquentielle / variance
        L = length(dN);
        
  
        n_win = 750;  
        n_ov  = 500;  % 75% overlap
        
        % NFFT : Zero-padding pour résolution spectrale fine
        % Résolution = Fs / NFFT = 5 / 2048 ≈ 0.0024 Hz
        current_nfft = TARGET_NFFT_PADDING;
        
        % Calculs spectraux
        [Pnn, f] = pwelch(dN, hamming(n_win), n_ov, current_nfft, FS_TARGET);
        [Pee, ~] = pwelch(dE, hamming(n_win), n_ov, current_nfft, FS_TARGET);
        [Pne, ~] = cpsd(dN, dE, hamming(n_win), n_ov, current_nfft, FS_TARGET);
        
        s1 = zeros(length(f), 1);
        
        % --- FILTRAGE PAR RATIO SVD ---
        for j = 1:length(f)
            Gyy = [Pnn(j), Pne(j); conj(Pne(j)), Pee(j)];
            [~, Sigma, ~] = svd(Gyy);
            
            val_s1 = Sigma(1,1);
            val_s2 = Sigma(2,2);
            
            if val_s2 > 0 && (val_s1 / val_s2) > DEFAULT_SVD_RATIO
                s1(j) = val_s1;
            else
                s1(j) = 1e-15;
            end
        end
        % ------------------------------
        
        svd_db = 10*log10(s1);
        
        %DETECTION
        mask_global = (f >= DEFAULT_FMIN) & (f <= DEFAULT_FMAX);
        if ~any(mask_global), continue; end
        
        sub_f_global = f(mask_global);
        sub_s_global = svd_db(mask_global);
        
        [max_val, idx_in_sub] = max(sub_s_global);
        f_at_max = sub_f_global(idx_in_sub);
        
        search_min = max(DEFAULT_FMIN, f_at_max - DEFAULT_DELTAF);
        search_max = min(DEFAULT_FMAX, f_at_max + DEFAULT_DELTAF);
        thresh_amp = max_val - DEFAULT_DYN;
        
        [pks, locs] = findpeaks(svd_db, f, 'MinPeakHeight', thresh_amp);
        mask_local = (locs >= search_min) & (locs <= search_max);
        final_locs = locs(mask_local);
        
        %STOCKAGE
        if ~isempty(final_locs)
            final_locs = sort(final_locs);
            nb = min(length(final_locs), NUM_MODES_MAX);
            all_peaks(k, 1:nb) = final_locs(1:nb);
        end
        
        %SNAPSHOT
        if k == idx_rand
            snap.valid = true;
            snap.date = dt;
            snap.f = f;
            snap.svd = svd_db;                                       
            snap.box = [search_min, search_max, thresh_amp, max_val]; 
            snap.pks = final_locs;                                   
            snap.amp = pks(mask_local);
            snap.n_win = n_win;
            snap.n_segments = floor((L - n_ov) / (n_win - n_ov));
        end
    end
    
    %% AFFICHAGE
    valid_idx = ~isnat(all_dates);
    d_valid = all_dates(valid_idx);
    peaks_valid = all_peaks(valid_idx, :);
    [d_valid, sort_idx] = sort(d_valid);
    peaks_valid = peaks_valid(sort_idx, :);
    
    if isempty(d_valid)
        warning("Pas de données valides pour %s", cfg.Name);
        continue;
    end
    
    figure('Name', sprintf("Suivi Cluster - %s", cfg.Name), 'Color', 'w', 'Position', [100+iTree*30, 100, 1100, 900]);
    tl = tiledlayout(3,1, 'Padding', 'compact');
    col_modes = lines(NUM_MODES_MAX);
    
    %BRUT
    nexttile; hold on; grid on;
    for m = 1:NUM_MODES_MAX
        scatter(d_valid, peaks_valid(:,m), 15, col_modes(m,:), 'filled', ...
            'DisplayName', sprintf('Pic %d', m));
    end
    ylabel('Fréquence (Hz)');
    title(sprintf('%s : Détection brute (Fmin=%.2f, Fmax=%.2f, Ratio=%.1f)', cfg.Name, DEFAULT_FMIN, DEFAULT_FMAX, DEFAULT_SVD_RATIO));
    ylim(DISPLAY_FREQ_LIMIT);
    legend('Location', 'eastoutside');
    
    %TENDANCE
    nexttile; hold on; grid on;
    for m = 1:NUM_MODES_MAX
        y_raw = peaks_valid(:,m);
        idx_ok = ~isnan(y_raw);
        if sum(idx_ok) < 10, continue; end
        
        t_sub = d_valid(idx_ok);
        y_sub = y_raw(idx_ok);
        y_smooth = smoothdata(y_sub, 'movmean', 10);
        
        t_plot = t_sub(1); y_plot = y_smooth(1);
        for j = 1:length(t_sub)-1
            if (t_sub(j+1) - t_sub(j)) > GAP_THRESHOLD
                t_plot = [t_plot; t_sub(j)+seconds(1); t_sub(j+1)]; %#ok<*AGROW>
                y_plot = [y_plot; NaN; y_smooth(j+1)];
            else
                t_plot = [t_plot; t_sub(j+1)];
                y_plot = [y_plot; y_smooth(j+1)];
            end
        end
        plot(t_plot, y_plot, 'LineWidth', 2, 'Color', col_modes(m,:), 'DisplayName', sprintf('Tendance Pic %d', m));
    end
    ylabel('Fréquence (Hz)');
    title('Tendance lissée du Cluster');
    ylim(DISPLAY_FREQ_LIMIT);
    
    %SNAPSHOT
    nexttile; 
    if snap.valid
        s_min = snap.box(1); s_max = snap.box(2);
        fill([s_min s_max s_max s_min], [-300 -300 300 300], ...
             [1 0.9 0.6], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
        hold on; grid on;
        plot(snap.f, snap.svd, 'k', 'LineWidth', 1.2);
        xline(s_min, '--', 'Color', [0.5 0.5 0.5]);
        xline(s_max, '--', 'Color', [0.5 0.5 0.5]);
        yline(snap.box(3), '--r', 'LineWidth', 1);
        if ~isempty(snap.pks)
            plot(snap.pks, snap.amp, 'rv', 'MarkerFaceColor','r');
            text(snap.pks, snap.amp + 1.5, string(round(snap.pks,3)), ...
                 'FontSize',8, 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
        end
        
        xlim(DISPLAY_FREQ_LIMIT);
        y_max_val = snap.box(4);
        ylim([y_max_val - Y_DYNAMIC_RANGE - 5, y_max_val + Y_MARGIN_UP]);
        
        title(sprintf('Snapshot : %s | Win=%d pts (%.1fs) | %d segments | Δf=%.4f Hz', ...
            string(snap.date), snap.n_win, snap.n_win/FS_TARGET, snap.n_segments, FS_TARGET/current_nfft)); 
        xlabel('Fréquence (Hz)'); ylabel('Amplitude (dB)');
    else
        text(0.5, 0.5, 'Pas de snapshot valide');
    end
end
disp("=== Terminé ===");

%% FONCTIONS LOCALES
function [d_out] = process_signal(raw, fs_in, fs_out, do_resample, order, target_len)
    d = double(raw(:));
    d = d - mean(d);
    
    current_len = length(d);
    if current_len < target_len
        d = [d; zeros(target_len - current_len, 1)];
    elseif current_len > target_len
        d = d(1:target_len);
    end
    
    % 1. Resampling
    if do_resample && (fs_in ~= fs_out)
        d = resample(d, fs_out, fs_in); 
    end
    
    % 2. Filtrage
    [b, a] = butter(order, [0.05 1.2] / (fs_out/2), 'bandpass');
    d_out = filtfilt(b, a, d);
end

function dt = parse_filename_date(fname)
    dt = NaT;
    [~,name,~] = fileparts(fname);
    tok = regexp(name, '(\d{4})(\d{2})(\d{2})[_-]?(\d{2})(\d{2})(\d{2})', 'tokens');
    if ~isempty(tok)
        v = tok{1};
        dt = datetime(str2double(v{1}), str2double(v{2}), str2double(v{3}), ...
                      str2double(v{4}), str2double(v{5}), str2double(v{6}));
    end
end