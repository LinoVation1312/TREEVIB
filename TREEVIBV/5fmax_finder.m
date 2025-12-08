clearvars; close all; clc;

%% PARAMETRES
mat_folder = "H:\Home\Documents\aaaTREEVIB\downsampling_filtered";
fs_target = 25;
f_search = [0.12 0.6];
Np = 5;                % Nombre de pics max à chercher
min_dist = 0.02;       % Distance minimale entre deux pics (Hz) - Zone d'exclusion
amp_thresh_ratio = 0.05; % Seuil d'amplitude relatif (5% du max global)
t_start = 20;
t_end = 280;
nfft_factor = 1;

%% CONFIGURATION ARBRES/ORIENTATIONS
configs = {'Arbre1_Nord', 'Arbre1_Est', 'Arbre2_Nord', 'Arbre2_Est'};
titles = {'Arbre 1 - Nord', 'Arbre 1 - Est', 'Arbre 2 - Nord', 'Arbre 2 - Est'};

files = dir(fullfile(mat_folder, "*.mat"));
if isempty(files), error("Aucun .mat trouvé"); end
nFiles = numel(files);

idx_rand = randi(nFiles);
fprintf("Fichier sélectionné pour affichage : %s (Index %d)\n", files(idx_rand).name, idx_rand);

%% BOUCLE PAR CONFIGURATION
for cfg = 1:length(configs)
    canal = configs{cfg};
    all_dates = NaT(nFiles,1);
    all_peaks = nan(nFiles, Np);
    rand_data = struct();
    
    for k = 1:nFiles
        fname = fullfile(files(k).folder, files(k).name);
        if mod(k,20)==0, fprintf("[%d/%d] %s  %s\n", k, nFiles, canal, files(k).name); end
        
        try
            S = load(fname);
        catch
            continue;
        end
        
        % --- Gestion Date ---
        dt = NaT;
        if isfield(S,'date_heure')
            try, dt = datetime(S.date_heure); catch, end
        end
        if isnat(dt), dt = parse_filename_date(files(k).name); end
        all_dates(k) = dt;
        
        % --- Chargement Signal ---
        if ~isfield(S, canal), continue; end
        x = S.(canal)(:);
        
        fs0 = 1000;
        if isfield(S,'fs'), fs0 = double(S.fs); end
        
        % --- Pre-processing ---
        x = resample(x, fs_target, fs0);
        
        idx_start = max(1, round(t_start * fs_target) + 1);
        idx_end = min(length(x), round(t_end * fs_target));
        if idx_start >= idx_end, continue; end
        x = x(idx_start:idx_end);
        
        x = detrend(x,'linear');
        [b,a] = butter(4, 1/(fs_target/2), 'low');
        x = filtfilt(b,a,x);
        
        % --- Calcul PSD (Welch) ---
        Nx = length(x);
        Nw = min(2^nextpow2(Nx), Nx);
        nfft = Nw * nfft_factor;
        win = hamming(Nw);
        ov = floor(Nw * 0.75);
        
        if Nw > Nx, continue; end
        
        [Pxx,f] = pwelch(x, win, ov, nfft, fs_target);
        
        % --- Restriction à la bande de recherche ---
        m = f >= f_search(1) & f <= f_search(2);
        ff = f(m);
        PP = Pxx(m); % Ceci est notre copie de travail pour la recherche
        PP_original = PP; % On garde l'original pour récupérer les vraies amplitudes
        
        if numel(ff) < 5, continue; end
        
        % === NOUVELLE METHODE : Iterative Peak Picking ===
        current_fpeaks = [];
        current_amps = [];
        
        max_amplitude = max(PP);
        abs_threshold = max_amplitude * amp_thresh_ratio;
        
        for iPeak = 1:Np
            % 1. Trouver le max global dans la bande actuelle
            [valMax, idxMax] = max(PP);
            
            % 2. Vérifier si ce max est significatif (au-dessus du bruit)
            if valMax < abs_threshold
                break; % On arrête si on ne trouve que du bruit
            end
            
            f_val = ff(idxMax);
            
            % 3. Stocker le pic
            current_fpeaks = [current_fpeaks; f_val];
            current_amps = [current_amps; valMax];
            
            % 4. Masquage Spectral (Spectral Masking)
            % On met à zéro les valeurs autour du pic trouvé pour ne pas
            % redétecter le même lobe ou un pic trop proche à la prochaine itération.
            mask_idx = abs(ff - f_val) <= min_dist;
            PP(mask_idx) = 0; 
        end
        
        % Tri des pics trouvés par fréquence croissante pour cohérence
        if ~isempty(current_fpeaks)
            [current_fpeaks, sort_idx] = sort(current_fpeaks);
            current_amps = current_amps(sort_idx);
            
            % Stockage
            Nfound = length(current_fpeaks);
            all_peaks(k, 1:Nfound) = current_fpeaks;
        end
        % =================================================
        
        if k == idx_rand
            rand_data.name = files(k).name;
            rand_data.date = dt;
            rand_data.f_full = f;
            rand_data.Pxx_full = Pxx;
            rand_data.ff = ff;
            rand_data.PP = PP_original; % Pour le plot, on veut le spectre non masqué
            rand_data.pks = current_amps;
            rand_data.locs = current_fpeaks;
        end
    end
    
    %% FIGURES PAR CONFIGURATION
    valid = ~isnat(all_dates);
    d_valid = all_dates(valid);
    peaks_valid = all_peaks(valid, :);
    
    figure('Name', titles{cfg}, 'Color', 'w', ...
           'Position', [100+cfg*50, 100+cfg*30, 1200, 900]);
    tiledlayout(3,1);
    
    %% --- GRAPH 1 : Modes bruts ---
    nexttile;
    hold on; grid on;
    colors = lines(Np);
    for i = 1:Np
        valid_mode = ~isnan(peaks_valid(:,i));
        if any(valid_mode)
            plot(d_valid(valid_mode), peaks_valid(valid_mode,i), ...
                'o', 'MarkerSize', 3, ...
                'MarkerFaceColor', colors(i,:), ...
                'Color', colors(i,:), ...
                'DisplayName', sprintf('Mode %d', i));
        end
    end
    ylabel("Fréquence [Hz]");
    title(sprintf("%s – Modes détectés (Méthode Itérative)", titles{cfg}));
    ylim([0.1 0.5])
    legend('Location','bestoutside');
    
%% --- GRAPH 2 : Modes lissés (Trait continu avec coupures) ---
nexttile;
hold on; grid on;

gap_threshold = hours(12); 

for i = 1:Np
    valid_mode = ~isnan(peaks_valid(:,i));
    
    if sum(valid_mode) > 5
        
        t_mode = d_valid(valid_mode);
        y_mode = peaks_valid(valid_mode, i);
        
        % 2. Lissage
        win_smooth = min(20, floor(sum(valid_mode)/5));
        if win_smooth < 5, win_smooth = 5; end
        y_smoothed = smoothdata(y_mode, 'movmean', win_smooth);
        

        t_plot = t_mode(1);
        y_plot = y_smoothed(1);
        
        for j = 1:length(t_mode)-1
            dt = t_mode(j+1) - t_mode(j);
            if dt > gap_threshold
                % S'il y a un trou, on ajoute un point NaN intermédiaire
                t_plot = [t_plot; t_mode(j) + seconds(1); t_mode(j+1)];
                y_plot = [y_plot; NaN; y_smoothed(j+1)];
            else
                % Sinon on continue la ligne normalement
                t_plot = [t_plot; t_mode(j+1)];
                y_plot = [y_plot; y_smoothed(j+1)];
            end
        end
        
        plot(t_plot, y_plot, '-', 'LineWidth', 1.5, ...
             'Color', colors(i,:), ...
             'DisplayName', sprintf('Mode %d', i));

    end
end

ylabel("Fréquence [Hz]");
title("Tendance lissée");
ylim([0.1 0.5])
legend('Location','bestoutside');
    
%% --- GRAPH 3 : Spectre aléatoire (Liste des fréquences en encart) ---
if isfield(rand_data, 'f_full')
    nexttile;
    
    % Préparation des données
    m_plot = rand_data.f_full >= 0 & rand_data.f_full <= 1;
    f_plot = rand_data.f_full(m_plot);
    Pxx_plot = rand_data.Pxx_full(m_plot);
    PP_dB = 10*log10(Pxx_plot);
    max_dB = max(PP_dB);
    min_dB_visu = max_dB - 30; 
    
    % 1. Trace du spectre
    plot(f_plot, PP_dB, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Signal'); 
    hold on;
    
    % 2. Limites de recherche
    xline(f_search(1), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'DisplayName', 'Limites');
    xline(f_search(2), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'HandleVisibility', 'off');
      
    if ~isempty(rand_data.locs)
        pks_dB = 10*log10(rand_data.pks);
        
        % 3. Trace des marqueurs rouges sur les pics (sans texte au-dessus)
        plot(rand_data.locs, pks_dB, 'rv', ...
             'MarkerFaceColor', 'r', 'MarkerSize', 7, ...
             'DisplayName', 'Pics détectés');
         
        % 4. CRÉATION DE LA LISTE DES FRÉQUENCES (Boîte "Légende")
        % On construit un tableau de cellules contenant le texte
        txt_legend = {'\bf Fréquences détectées :'}; % \bf pour le gras
        for i = 1:length(rand_data.locs)
            txt_legend{end+1} = sprintf('Mode %d : %.3f Hz', i, rand_data.locs(i));
        end
        
        % Affichage de la boîte en haut à droite (coordonnées normalisées)
        text(0.98, 0.96, txt_legend, ...
             'Units', 'normalized', ...       % Position relative à la taille de la figure (0 à 1)
             'VerticalAlignment', 'top', ...  % Ancrage haut
             'HorizontalAlignment', 'right',...% Ancrage droite
             'FontSize', 9, ...
             'BackgroundColor', [1 1 1 0.9], ... % Fond blanc légèrement transparent
             'EdgeColor', 'k', ...            % Bordure noire
             'Margin', 5);                    % Marge intérieure
    end
    
    grid on;
    xlabel('Fréquence [Hz]');
    ylabel('DSP [dB]');
    title(sprintf('Spectre – %s', string(rand_data.date)));
    xlim([0 1]);
    
    if ~isempty(max_dB)
        ylim([min_dB_visu, max_dB + 5]); % Marge réduite car plus de texte au-dessus
    end
    
    % Légende standard pour les courbes uniquement, placée en bas à gauche
    legend('Location','southwest');
end
end
disp("=== FIN ===");

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