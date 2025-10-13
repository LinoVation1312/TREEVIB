% Lecture et traitement de fichiers TDMS avec filtrage
clear all; close all; clc;

%% ========== PARAMÈTRES CONFIGURABLES ==========
% Paramètres du fichier
tdms_path ="V:\MONITORING_ARBRES\PEUPLIERS a partir du 03 mars (a mettre a jour chaque mois)\mesure_2025_05_28_23_13_59.tdms";

% Paramètres du signal
fs_original = 1000; 
fs = 40;             

% Paramètres de filtrage
f_dc = 0.1;         % HP [Hz]
fc = 0.5;            % LP [Hz]
hp_taps = 2001;      % taps FIR HP
lp_order = 2;        % Ordre du filtre Butterworth

% Paramètres de pwelch
pwelch_window = 300*fs;
pwelch_overlap = 0 * pwelch_window;
pwelch_nfft = 2^nextpow2(pwelch_window);
pwelch_freq_range = [0 1]; 

% Organisation par arbre
trees = struct(...
    'Arbre1', {{'Arbre1_Nord', 'Arbre1_Est'}}, ...
    'Arbre2', {{'Arbre2_Nord', 'Arbre2_Est'}} ...
);

% ================================================

%% Vérification de l'existence du fichier
if ~isfile(tdms_path)
    error("Le fichier %s est introuvable.", tdms_path);
end

%% Lecture des métadonnées TDMS
fprintf('\n=== Informations TDMS ===\n');
info = tdmsinfo(tdms_path);

%% Sélection d'un groupe TDMS (via ChannelList)
if ~isprop(info, 'ChannelList')
    error('tdmsinfo ne fournit pas ChannelList. Impossible de lister les groupes.');
end

chanlist = info.ChannelList; % table avec colonnes: ChannelGroupName, ChannelName, etc.
if ~all(ismember({'ChannelGroupName','ChannelName'}, chanlist.Properties.VariableNames))
    error('ChannelList ne contient pas ChannelGroupName/ChannelName dans cette version.');
end

% Liste des groupes disponibles
groupes = unique(chanlist.ChannelGroupName);
fprintf('\n=== Groupes TDMS disponibles ===\n');
disp(groupes);

% Sélection du groupe par interface ou console
try
    [idx, ok] = listdlg('PromptString','Sélectionner un groupe:', ...
                        'SelectionMode','single', ...
                        'ListString', cellstr(groupes), ...
                        'Name','Sélection du groupe');
    if ~ok, error('Sélection annulée.'); end
    selected_group = groupes(idx);
catch
    fprintf('Interface indisponible. Sélection via console.\n');
    for gi = 1:numel(groupes)
        fprintf('  [%d] %s\n', gi, groupes{gi});
    end
    gi = input('Entrez l''indice du groupe à analyser: ');
    if ~(isnumeric(gi) && gi>=1 && gi<=numel(groupes))
        error('Indice invalide.');
    end
    selected_group = groupes(gi);
end

fprintf('\nGroupe sélectionné: %s\n', selected_group);

%% Canaux appartenant au groupe sélectionné
in_group = strcmp(chanlist.ChannelGroupName, selected_group);
channels_in_group = chanlist.ChannelName(in_group);

fprintf('Canaux dans le groupe sélectionné:\n');
disp(channels_in_group);

%% Lecture uniquement du groupe choisi
data_table = tdmsread(tdms_path, ...
                      ChannelGroupName = selected_group, ...
                      ChannelNames = channels_in_group);

% Déballer si c'est une cellule contenant une seule table
if iscell(data_table) && numel(data_table) == 1 && istable(data_table{1})
    data_table = data_table{1};
end

fprintf('\n=== Données chargées pour le groupe %s ===\n', selected_group);

if istable(data_table)
    disp(data_table.Properties.VariableNames);
elseif iscell(data_table)
    fprintf('Résultat = cellule (%d éléments)\n', numel(data_table));
    for i = 1:numel(data_table)
        if istable(data_table{i})
            fprintf(' -> Table à l''index %d avec colonnes:\n', i);
            disp(data_table{i}.Properties.VariableNames);
        else
            fprintf(' -> Élément %d de type %s\n', i, class(data_table{i}));
        end
    end
else
    fprintf('Type inattendu: %s\n', class(data_table));
    disp(data_table);
end

%% Séparer canaux temporels et DSP
is_dsp = endsWith(channels_in_group, "_DSP");
channels_time = channels_in_group(~is_dsp);
channels_dsp  = channels_in_group(is_dsp);

if isempty(channels_time)
    error('Aucun canal temporel présent dans la table pour ce groupe.');
end


%% Affichage des paramètres
fprintf('\n=== Paramètres de traitement ===\n');
fprintf('Fréquence originale : %d Hz\n', fs_original);
fprintf('Fréquence après resampling : %d Hz\n', fs);
fprintf('Taille de fenêtre pwelch : %d\n', pwelch_window);
fprintf('Recouvrement : %d échantillons\n', pwelch_overlap);
fprintf('NFFT : %d\n', pwelch_nfft);
fprintf('Résolution fréquentielle : %.4f Hz\n', fs/pwelch_nfft);

%% Adapter la structure trees au groupe sélectionné
tree_names = fieldnames(trees);
for ti = 1:numel(tree_names)
    chs = trees.(tree_names{ti});
    trees.(tree_names{ti}) = intersect(chs, channels_time, 'stable');
end

% Supprimer les arbres vides
tree_names = fieldnames(trees);
empty_mask = cellfun(@(t) isempty(trees.(t)), tree_names);
trees = rmfield(trees, tree_names(empty_mask));
tree_names = fieldnames(trees);

if isempty(tree_names)
    % Si aucun arbre valide, créer un arbre générique avec les canaux du groupe choisi
    trees = struct('Groupe', cellstr(channels_time(:)'));
    tree_names = fieldnames(trees);
    fprintf('Aucun arbre prédéfini valide pour ce groupe. Analyse générique par canaux.\n');
end

colors = [0 0.447 0.741; 0.85 0.33 0.1];  % Bleu pour Nord, Orange pour Est

%% Traitement par arbre
output_root = fullfile(pwd);

for tree_idx = 1:length(tree_names)
    tree_name = tree_names{tree_idx};
    channels = trees.(tree_name);

    fprintf('\n========================================\n');
    fprintf('Traitement de : %s (groupe %s)\n', tree_name, selected_group);
    fprintf('========================================\n');

    % Préparation du stockage des données
    data_processed = struct();

    % Traitement de chaque canal
    for dir_idx = 1:length(channels)
        channel_name = channels{dir_idx};

        if ~ismember(channel_name, data_table.Properties.VariableNames)
            fprintf('ATTENTION : Canal %s non trouvé, ignoré.\n', channel_name);
            continue;
        end

        fprintf('\n  -> %s\n', channel_name);

        % Extraction des données à fs_original
        data_1000Hz = data_table.(channel_name);
        data_1000Hz = double(data_1000Hz(:));
        N_original = length(data_1000Hz);

        % Resampling à fs Hz
        N_target = round(N_original * fs / fs_original);
        data = resample(data_1000Hz, N_target, N_original);
        clear data_1000Hz;

        fprintf('     Points à %d Hz : %d\n', fs, length(data));
        fprintf('     Durée : %.2f s\n', length(data) / fs);

        % Calcul du temps
        N = length(data);
        t = (0:N-1) / fs;

        % PSD brute
        [Pxx, f_psd] = pwelch(data, pwelch_window, pwelch_overlap, pwelch_nfft, fs);

        % Filtrage
        data_bp = filtfilt_cascade(data, fs, f_dc, fc, hp_taps, lp_order);

        % PSD filtrée
        [Pxx_bp, f_psd_bp] = pwelch(data_bp, pwelch_window, pwelch_overlap, pwelch_nfft, fs);

        % Stockage
        data_processed.(channel_name) = struct(...
            't', t, ...
            'data', data, ...
            'data_bp', data_bp, ...
            'Pxx', Pxx, ...
            'f_psd', f_psd, ...
            'Pxx_bp', Pxx_bp, ...
            'f_psd_bp', f_psd_bp ...
        );
    end



    %% Figure
    figure('Position', [100, 100, 1400, 900]);
    figure_title = sprintf('Analyse %s (groupe %s) - fs = %d Hz', tree_name, selected_group, fs);
    sgtitle(figure_title, 'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');

    % Tracé 1: signal brut
    subplot(2,2,1); hold on; grid on;
    for dir_idx = 1:length(channels)
        channel_name = channels{dir_idx};
        if ~isfield(data_processed, channel_name), continue; end
        d = data_processed.(channel_name);
        c = colors(1 + mod(dir_idx-1, size(colors,1)), :);
        plot(d.t, d.data, 'Color', c, 'DisplayName', channel_name);
    end
    xlabel('Temps [s]'); ylabel('Amplitude'); title('Signal brut'); legend('Interpreter','none');

    % Tracé 2: signal filtré
    subplot(2,2,2); hold on; grid on;
    for dir_idx = 1:length(channels)
        channel_name = channels{dir_idx};
        if ~isfield(data_processed, channel_name), continue; end
        d = data_processed.(channel_name);
        c = colors(1 + mod(dir_idx-1, size(colors,1)), :);
        plot(d.t, d.data_bp, 'Color', c, 'DisplayName', channel_name);
    end
    xlabel('Temps [s]'); ylabel('Amplitude'); title('Signal filtré (HP+LP)'); legend('Interpreter','none');

    % Tracé 3: PSD brute
    subplot(2,2,3); hold on; grid on;
    for dir_idx = 1:length(channels)
        channel_name = channels{dir_idx};
        if ~isfield(data_processed, channel_name), continue; end
        d = data_processed.(channel_name);
        c = colors(1 + mod(dir_idx-1, size(colors,1)), :);
        plot(d.f_psd, d.Pxx, 'Color', c, 'DisplayName', channel_name);xlim(pwelch_freq_range);
    end
     xlabel('Fréquence [Hz]'); ylabel('PSD [m/s²/Hz]'); title('PSD brute'); legend('Interpreter','none');
    

    % Tracé 4: PSD filtrée
    subplot(2,2,4); hold on; grid on;
    for dir_idx = 1:length(channels)
        channel_name = channels{dir_idx};
        if ~isfield(data_processed, channel_name), continue; end
        d = data_processed.(channel_name);
        c = colors(1 + mod(dir_idx-1, size(colors,1)), :);
        plot(d.f_psd_bp, d.Pxx_bp, 'Color', c, 'DisplayName', channel_name);xlim(pwelch_freq_range);
    end
     xlabel('Fréquence [Hz]'); ylabel('[m/s²/Hz]'); title('PSD filtrée'); legend('Interpreter','none');
    
end

fprintf('\n=== Terminé ===\n');

%% ========= Fonctions auxiliaires =========
function y = filtfilt_cascade(x, fs, f_dc, fc, hp_taps, lp_order)
% Filtre passe-haut FIR (soustraction DC) + passe-bas Butterworth, appliqués avec filtfilt

    % HP FIR (type I, atténuation DC)
    nyq = fs/2;
    if f_dc <= 0
        y_hp = x;
    else
        Wn_hp = f_dc/nyq;
        % Conception FIR HP via design de bandstop sur DC (using fir1 with highpass)
        b_hp = fir1(hp_taps-1, Wn_hp, 'high', hamming(hp_taps), 'scale');
        y_hp = filtfilt(b_hp, 1, x);
    end

    % LP Butterworth
    if fc >= nyq
        y = y_hp;
        return;
    end
    Wn_lp = fc/nyq;
    [b_lp, a_lp] = butter(lp_order, Wn_lp, 'low');
    y = filtfilt(b_lp, a_lp, y_hp);
end
