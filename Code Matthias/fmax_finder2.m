%% ============================================================
%  analyse_cluster_3_plots.m
%  Génère 3 graphiques distincts alignés verticalement :
%  1. Fréquence vs Température
%  2. Fréquence vs Pression
%  3. Fréquence vs Vent
% ============================================================
clearvars; close all; clc;

%% 1. PARAMÈTRES
mat_folder = "H:\Home\Documents\1-5A_ProjetTreeVib\data de Mars à Novembre 2025 - fs=25Hz";
fs_target = 25;                 
f_search = [0.05 0.6];           

% --- Paramètres Arbres ---
t_start = 50; t_end = 250;                    
min_window = 512; max_window = 2048;
smoothing_freq = 50;  

% --- Paramètres Météo ---
% Fenêtre de lissage (jours)
window_meteo = 30; 


%% 2. TRAITEMENT ARBRES
files = dir(fullfile(mat_folder, "*.mat"));
if isempty(files)
    error("Aucun .mat trouvé dans %s", mat_folder);
end
fprintf("Fichiers trouvés : %d\n", numel(files));

channel_names = {'Arbre1_Nord', 'Arbre1_Est', 'Arbre2_Nord', 'Arbre2_Est'};
line_colors   = {'r', 'b', 'm', 'g'}; 

all_dates = NaT(numel(files),1);
for i = 1:length(channel_names)
    results.(channel_names{i}) = nan(numel(files), 1);
end

for k = 1:numel(files)
    fname = fullfile(files(k).folder, files(k).name);
    try S = load(fname); catch, continue; end
    
    dt = NaT;
    if isfield(S,'date_heure'), try dt = datetime(S.date_heure); catch, end; end
    if isnat(dt), dt = parse_filename_date(files(k).name); end
    all_dates(k) = dt;
    
    if isfield(S,'fs'), fs0 = double(S.fs); else, fs0 = 1000; end
    fprintf("[%d/%d] %s ... ", k, numel(files), files(k).name);

    for c = 1:length(channel_names)
        chan = channel_names{c};
        if isfield(S, chan)
            x_raw = S.(chan)(:);
            if fs0 ~= fs_target, x = resample(x_raw, fs_target, fs0); else, x = x_raw; end
            idx_start = round(t_start * fs_target) + 1; idx_end = round(t_end * fs_target);
            if idx_start < idx_end && idx_start < length(x)
                x = x(max(1,idx_start):min(length(x),idx_end));
                x = detrend(x,'linear');
                [b,a] = butter(4, 1/(fs_target/2), 'low'); x = filtfilt(b,a,x);
                Nw = min(max_window, max(min_window, floor(length(x)/2)));
                [Pxx, f] = pwelch(x, hamming(Nw), floor(Nw/2), 2^nextpow2(Nw), fs_target);
                mask = f >= f_search(1) & f <= f_search(2);
                ff = f(mask); PP = Pxx(mask);
                if sum(PP) > 0 && numel(ff) > 5
                    results.(chan)(k) = sum(ff .* PP) / sum(PP);
                end
            end
        end
    end
    fprintf("OK\n");
end
valid = ~isnat(all_dates);
dates_valid = all_dates(valid);

%% 3. MÉTÉO : Chargement TEMP, VENT, PRESSION
disp('--- Chargement Meteostat ---');
[fileMeteo, pathMeteo] = uigetfile({'*.csv;*.xlsx', 'Fichiers Météo'}, 'Sélectionner le fichier Meteostat');
TT_Meteo_Plot = table(); 

if ~isequal(fileMeteo, 0)
    fullPathMeteo = fullfile(pathMeteo, fileMeteo);
    opts = detectImportOptions(fullPathMeteo);
    opts.VariableNamingRule = 'preserve';
    tableMeteo = readtable(fullPathMeteo, opts);
    
    % Extraction colonnes (selon votre capture)
    d_raw = tableMeteo{:, 1}; % Date
    t_raw = tableMeteo{:, 2}; % Température
    v_raw = tableMeteo{:, 3}; % Vent
    p_raw = tableMeteo{:, 4}; % Pression
    
    try d_raw = datetime(d_raw, 'InputFormat', 'yyyy-MM-dd HH:mm:ss'); catch, try d_raw = datetime(d_raw); catch, end, end
    
    TT_MeteoStat = timetable(d_raw, t_raw, v_raw, p_raw, 'VariableNames', {'Temp', 'Vent', 'Pression'});
    
    % Moyenne journalière + Lissage Gaussien
    TT_Meteo_Plot = retime(TT_MeteoStat, 'daily', 'mean');
    
    if height(TT_Meteo_Plot) > 0
        TT_Meteo_Plot.Temp     = smoothdata(TT_Meteo_Plot.Temp, 'gaussian', window_meteo);
        TT_Meteo_Plot.Vent     = smoothdata(TT_Meteo_Plot.Vent, 'gaussian', window_meteo);
        TT_Meteo_Plot.Pression = smoothdata(TT_Meteo_Plot.Pression, 'gaussian', window_meteo);
        disp(['Météo chargée et lissée (Fenêtre ' num2str(window_meteo) ' jours).']);
    end
else
    disp('Pas de fichier météo.');
end

%% 4. FIGURE : 3 TILES (Temp, Pression, Vent)
figure('Name', 'Analyses Météo Distinctes', 'Color', 'w', 'Position', [50 50 1000 900]);
tlo = tiledlayout(3,1, 'TileSpacing', 'compact'); 

% --- 1. FREQ vs TEMP ---
ax1 = nexttile;
yyaxis left; hold on;
for c=1:4, plot(dates_valid, smoothdata(results.(channel_names{c})(valid),'movmean',smoothing_freq), 'Color', line_colors{c}, 'LineWidth',1.2); end
ylabel('Fréquence [Hz]'); ylim([0.15 0.4]); set(gca,'YColor','k');
yyaxis right; 
if ~isempty(TT_Meteo_Plot)
    plot(TT_Meteo_Plot.d_raw, TT_Meteo_Plot.Temp, 'Color', [0.9 0.4 0.1], 'LineWidth', 2);
    ylabel('Température (°C)'); ylim([0 35]); set(gca,'YColor',[0.9 0.4 0.1]);
end
title('1. Fréquence vs Température'); grid on;

% --- 2. FREQ vs PRESSION ---
ax2 = nexttile;
yyaxis left; hold on;
for c=1:4, plot(dates_valid, smoothdata(results.(channel_names{c})(valid),'movmean',smoothing_freq), 'Color', line_colors{c}, 'LineWidth',1.2); end
ylabel('Fréquence [Hz]'); ylim([0.15 0.4]); set(gca,'YColor','k');
yyaxis right; 
if ~isempty(TT_Meteo_Plot)
    plot(TT_Meteo_Plot.d_raw, TT_Meteo_Plot.Pression, 'Color', [0.5 0 0.5], 'LineWidth', 2);
    ylabel('Pression (hPa)'); 
    % Zoom automatique centré sur la pression moyenne
    if ~all(isnan(TT_Meteo_Plot.Pression))
        p_mean = mean(TT_Meteo_Plot.Pression, 'omitnan');
        ylim([p_mean-15, p_mean+15]); 
    end
    set(gca,'YColor',[0.5 0 0.5]);
end
title('2. Fréquence vs Pression Atmosphérique'); grid on;

% --- 3. FREQ vs VENT ---
ax3 = nexttile;
yyaxis left; hold on;
for c=1:4, plot(dates_valid, smoothdata(results.(channel_names{c})(valid),'movmean',smoothing_freq), 'Color', line_colors{c}, 'LineWidth',1.2); end
ylabel('Fréquence [Hz]'); ylim([0.15 0.4]); set(gca,'YColor','k');
yyaxis right; 
if ~isempty(TT_Meteo_Plot)
    plot(TT_Meteo_Plot.d_raw, TT_Meteo_Plot.Vent, 'Color', [0 0.6 0.7], 'LineWidth', 2);
    ylabel('Vent (km/h)'); 
    ylim([0 max(TT_Meteo_Plot.Vent)*1.2]); % 0 à Max
    set(gca,'YColor',[0 0.6 0.7]);
end
title('3. Fréquence vs Vent Moyen'); grid on;
xlabel('Date');

% Liaison Axes X
linkaxes([ax1, ax2, ax3], 'x');
datetick(ax3, 'x', 'dd-mmm', 'keeplimits');
if ~isempty(TT_Meteo_Plot), xlim([min(dates_valid) max(TT_Meteo_Plot.d_raw)]); end

% Légende unique en haut
legend({'A1 Nord', 'A1 Est', 'A2 Nord', 'A2 Est', 'Donnée Météo'}, 'Location', 'northoutside', 'Orientation', 'horizontal');

%% Fonction Date
function dt = parse_filename_date(fname)
    dt = NaT; [~,name,~] = fileparts(fname);
    tok = regexp(name, '(\d{4})(\d{2})(\d{2})[_-]?(\d{2})(\d{2})(\d{2})', 'tokens');
    if ~isempty(tok), v = tok{1}; dt = datetime(str2double(v{1}),str2double(v{2}),str2double(v{3}),str2double(v{4}),str2double(v{5}),str2double(v{6})); end
end