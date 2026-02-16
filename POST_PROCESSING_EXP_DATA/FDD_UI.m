classdef FDD_UI < handle
    properties (Access = public)
        %% ========== CONFIGURATION  ==========
        MAT_FOLDER = "H:\Home\Documents\aaaTREEVIB\new_5Hz";
        FS_TARGET = 5;              
        DEFAULT_FMIN = 0.1;          
        DEFAULT_FMAX = 0.75;
        
        % --- PARAMETRES RESTAURÉS & MODIFIÉS ---
        DEFAULT_DELTAF = 0.1;        % Largeur fenêtre (Hz) autour du pic
        DEFAULT_DYN = 10;            % Seuil relatif (dB)
        
        % MODIFICATION ICI : Ratio physique (10 = 10dB) au lieu de 1e+16
        DEFAULT_SVD_RATIO = 4;     
        
        DEFAULT_NFFT = 1500;         % Points pour la finesse du tracé (Zero-padding)
        BUTTER_ORDER = 2;            
        DISPLAY_FREQ_LIMIT = [0 1];  
        
        % Paramètres de dynamique d'affichage (Y)
        Y_MARGIN_UP = 5;             
        Y_DYNAMIC_RANGE = 30;        
        
        % Couleurs
        ColNord = [0 0.4470 0.7410]; 
        ColEst  = [0.8500 0.3250 0.0980]; 
        ColSVD  = [0.15 0.15 0.15];  
        
        TreesConfig = struct('Arbre1', {{'Arbre1_Nord', 'Arbre1_Est'}}, ...
                             'Arbre2', {{'Arbre2_Nord', 'Arbre2_Est'}})
    end
    
    properties (Access = private)
        UIFigure, FileListBox
        AxTemp1, AxTemp2
        AxFreq1, AxFreq2
        FMinSpin, FMaxSpin, NFFTEdit
        % Contrôles restaurés :
        DeltaFSpin, DynEdit, RatioSpin 
        DataStruct
    end
    
    methods
        function obj = FDD_UI()
            clc; close all; 
            obj.createUI();
            obj.loadFileList();
        end
        
        function createUI(obj)
            obj.UIFigure = uifigure('Name', 'Analyse FDD : Ratio + Fenêtrage Dynamique', 'Position', [100 100 1300 850]);
            gl = uigridlayout(obj.UIFigure, [1 2]);
            gl.ColumnWidth = {260, '1x'};
            
            p = uipanel(gl, 'Title', 'Configuration');
            uilabel(p, 'Position', [10 650 200 22], 'Text', 'Mesures disponibles :');
            
            obj.FileListBox = uilistbox(p, 'Position', [10 450 230 180], ...
                'ValueChangedFcn', @(src,event) obj.loadFile());
            
            y = 410;
            uilabel(p, 'Position', [10 y 120 22], 'Text', 'Plage fmin (Hz)');
            obj.FMinSpin = uieditfield(p, 'numeric', 'Position', [140 y 60 22], 'Value', obj.DEFAULT_FMIN);
            
            y = y-30;
            uilabel(p, 'Position', [10 y 120 22], 'Text', 'Plage fmax (Hz)');
            obj.FMaxSpin = uieditfield(p, 'numeric', 'Position', [140 y 60 22], 'Value', obj.DEFAULT_FMAX);
            
            y = y-30;
            uilabel(p, 'Position', [10 y 120 22], 'Text', 'Ratio S1/S2 Min');
            obj.RatioSpin = uieditfield(p, 'numeric', 'Position', [140 y 60 22], 'Value', obj.DEFAULT_SVD_RATIO);
            
            % --- UI DELTA F & DYN ---
            y = y-30;
            uilabel(p, 'Position', [10 y 120 22], 'Text', 'Delta f (+/- Hz)');
            obj.DeltaFSpin = uieditfield(p, 'numeric', 'Position', [140 y 60 22], 'Value', obj.DEFAULT_DELTAF);
            
            y = y-30;
            uilabel(p, 'Position', [10 y 120 22], 'Text', 'Dyn. (dB)');
            obj.DynEdit = uieditfield(p, 'numeric', 'Position', [140 y 60 22], 'Value', obj.DEFAULT_DYN);
            
            y = y-30;
            uilabel(p, 'Position', [10 y 120 22], 'Text', 'NFFT (Padding)');
            obj.NFFTEdit = uieditfield(p, 'numeric', 'Position', [140 y 60 22], 'Value', obj.DEFAULT_NFFT);
            
            uibutton(p, 'Position', [10 180 230 40], 'Text', 'CALCULER FDD', 'BackgroundColor', [.1 .4 .8], 'FontColor', 'w', 'FontWeight', 'bold', 'ButtonPushedFcn', @(src,event) obj.runFDD());
            uibutton(p, 'Position', [10 130 230 30], 'Text', 'RESET INTERFACE', 'ButtonPushedFcn', @(src,event) obj.resetApp());
            
            gr = uigridlayout(gl, [2 2]);
            obj.AxTemp1 = uiaxes(gr); obj.AxTemp2 = uiaxes(gr);
            obj.AxFreq1 = uiaxes(gr); obj.AxFreq2 = uiaxes(gr);
            
            obj.applyTitles();
            grid([obj.AxTemp1, obj.AxTemp2, obj.AxFreq1, obj.AxFreq2], 'on');
        end
        
        function applyTitles(obj)
            title(obj.AxTemp1, 'Arbre 1 : Temporel'); title(obj.AxTemp2, 'Arbre 2 : Temporel');
            title(obj.AxFreq1, 'Arbre 1 : SVD & Détection'); title(obj.AxFreq2, 'Arbre 2 : SVD & Détection');
            ylabel(obj.AxFreq1, 'dB/Hz'); ylabel(obj.AxFreq2, 'dB/Hz');
            xlabel(obj.AxFreq1, 'Fréquence (Hz)'); xlabel(obj.AxFreq2, 'Fréquence (Hz)');
        end
        
        function loadFileList(obj)
            if ~exist(obj.MAT_FOLDER, 'dir')
                obj.FileListBox.Items = "Dossier introuvable"; return; 
            end
            files = dir(fullfile(obj.MAT_FOLDER, '*.mat'));
            if isempty(files)
                obj.FileListBox.Items = "Aucun fichier .mat"; return;
            end
            [~, idx] = sort({files.name}); files = files(idx);
            rawNames = string({files.name});
            displayNames = rawNames; 
            for i = 1:length(rawNames)
                tokens = regexp(rawNames(i), '(\d{8})_(\d{6})', 'tokens', 'once');
                if ~isempty(tokens)
                    try
                        dt = datetime(string(tokens{1}) + string(tokens{2}), 'InputFormat', 'yyyyMMddHHmmss');
                        displayNames(i) = string(dateshift(dt, 'start', 'minute', 'nearest'), 'dd/MM - HH''h''mm');
                    catch; end
                end
            end
            obj.FileListBox.Items = displayNames;      
            obj.FileListBox.ItemsData = rawNames;      
        end
        
        function loadFile(obj)
            selected = obj.FileListBox.Value;
            if isempty(selected) || strcmp(selected, 'Dossier introuvable'), return; end
            try
                obj.DataStruct = load(fullfile(obj.MAT_FOLDER, selected));
                fprintf('Chargement : %s\n', selected);
                obj.plotPreview();
            catch ME
                uialert(obj.UIFigure, ME.message, 'Erreur');
            end
        end
        
        function plotPreview(obj)
            obj.cleanAxes(obj.AxTemp1); obj.cleanAxes(obj.AxTemp2);
            obj.cleanAxes(obj.AxFreq1); obj.cleanAxes(obj.AxFreq2);
            obj.applyTitles();
            obj.plotTreeData(obj.AxTemp1, obj.AxFreq1, 'Arbre1');
            obj.plotTreeData(obj.AxTemp2, obj.AxFreq2, 'Arbre2');
        end
        
        function plotTreeData(obj, axT, axF, treeKey)
            hold(axT, 'on'); hold(axF, 'on');
            names = obj.TreesConfig.(treeKey);
            if ~isfield(obj.DataStruct, names{1}) || ~isfield(obj.DataStruct, names{2})
                text(axT, 0.5, 0.5, 'Données manquantes', 'HorizontalAlignment', 'center'); return;
            end
            
            [t, dN] = obj.process_signal(obj.DataStruct.(names{1}));
            [~, dE] = obj.process_signal(obj.DataStruct.(names{2}));
            
            plot(axT, t, dN, 'Color', obj.ColNord, 'DisplayName', 'Nord');
            plot(axT, t, dE, 'Color', obj.ColEst, 'DisplayName', 'Est');
            legend(axT, 'show', 'FontSize', 7);
            
            nfft = obj.NFFTEdit.Value;
            % Prévisualisation standard (Welch classique)
            [Pnn, f] = pwelch(dN, hamming(nfft), nfft/2, nfft, obj.FS_TARGET);
            [Pee, ~] = pwelch(dE, hamming(nfft), nfft/2, nfft, obj.FS_TARGET);
            
            plot(axF, f, 10*log10(Pnn), 'Color', [obj.ColNord 0.4], 'LineWidth', 0.8, 'DisplayName', 'DSP Nord');
            plot(axF, f, 10*log10(Pee), 'Color', [obj.ColEst 0.4], 'LineWidth', 0.8, 'DisplayName', 'DSP Est');
            
            xlim(axF, obj.DISPLAY_FREQ_LIMIT);
            ylim(axF); 
        end
        
        function runFDD(obj)
            if isempty(obj.DataStruct), return; end
            obj.computeSVD(obj.AxFreq1, 'Arbre1');
            obj.computeSVD(obj.AxFreq2, 'Arbre2');
        end
        
        function computeSVD(obj, ax, treeKey)
            % Nettoyage
            delete(findall(ax, 'Tag', 'SVD_Plot'));
            delete(findall(ax, 'Type', 'text'));
            delete(findall(ax, 'Type', 'fill')); 
            delete(findall(ax, 'Tag', 'PeakMarker'));
            
            names = obj.TreesConfig.(treeKey);
            if ~isfield(obj.DataStruct, names{1}), return; end
            [~, dN] = obj.process_signal(obj.DataStruct.(names{1}));
            [~, dE] = obj.process_signal(obj.DataStruct.(names{2}));
            
            % Sécurité dimensions
            nMin = min(length(dN), length(dE));
            dN = dN(1:nMin); dE = dE(1:nMin);
            
            % Paramètres UI
            nfft_user = obj.NFFTEdit.Value;
            ratio_min = obj.RatioSpin.Value;
            f_min = obj.FMinSpin.Value; 
            f_max = obj.FMaxSpin.Value;
            delta_f = obj.DeltaFSpin.Value;
            dyn_val = obj.DynEdit.Value;
            
            % --- CALCUL FORCE POUR 3 MOYENNES (Welch 50% Overlap) ---
            % L = Win + (2 * 0.5 * Win) = 2 * Win  => Win = L / 2
            L = length(dN);
            n_win = floor(L / 2);         
            n_overlap = floor(n_win / 2); 
            
            % NFFT doit être >= fenêtre. 
            % On prend le max entre le NFFT désiré (finesse) et la fenêtre requise (physique)
            nfft = max(nfft_user, n_win);
            
            [Pnn, f] = pwelch(dN, hamming(n_win), n_overlap, nfft, obj.FS_TARGET);
            [Pee, ~] = pwelch(dE, hamming(n_win), n_overlap, nfft, obj.FS_TARGET);
            [Pne, ~] = cpsd(dN, dE, hamming(n_win), n_overlap, nfft, obj.FS_TARGET);
            
            s1 = zeros(length(f), 1);
            
            % 1. BOUCLE SVD + FILTRE RATIO
            for j = 1:length(f)
                Gyy = [Pnn(j), Pne(j); conj(Pne(j)), Pee(j)];
                [~, S, ~] = svd(Gyy); 
                val_s1 = S(1,1); val_s2 = S(2,2);
                
                % Le coeur du filtrage : On vérifie le ratio physique
                if val_s2 > 0 && (val_s1 / val_s2) > ratio_min
                    s1(j) = val_s1;
                else
                    s1(j) = 1e-15; % Bruit de fond écrasé
                end
            end
            
            psd_db = 10*log10(s1);
            plot(ax, f, psd_db, 'Color', obj.ColSVD, 'LineWidth', 2, 'DisplayName', 'SVD Filtrée', 'Tag', 'SVD_Plot');
            
            % 2. FENETRAGE DYNAMIQUE + SEUILLAGE
            in_range = (f >= f_min) & (f <= f_max);
            if any(in_range)
                sub_psd = psd_db(in_range);
                sub_f = f(in_range);
                
                % A. Trouver le maximum dans la zone globale
                [max_val, idx_in_sub] = max(sub_psd);
                f_at_max = sub_f(idx_in_sub);
                
                % B. Définir la fenêtre jaune autour du max
                search_min = max(f_min, f_at_max - delta_f);
                search_max = min(f_max, f_at_max + delta_f);
                
                % C. Définir le seuil d'amplitude dynamique
                thresh_abs = max_val - dyn_val;
                
                % D. Détection des pics UNIQUEMENT dans cette fenêtre & au-dessus du seuil
                [pks, locs] = findpeaks(psd_db, f, 'MinPeakHeight', thresh_abs);
                
                % Filtre final : garder uniquement les pics DANS la fenêtre
                mask_window = (locs >= search_min) & (locs <= search_max);
                final_locs = locs(mask_window);
                final_pks  = pks(mask_window);
                
                % E. Affichage (Pics + Zone jaune)
                y_max_view = max_val + obj.Y_MARGIN_UP;
                y_min_view = y_max_view - obj.Y_DYNAMIC_RANGE;
                ylim(ax, [y_min_view, y_max_view]);
                
                % Dessin zone jaune
                fill(ax, [search_min search_max search_max search_min], ...
                     [y_min_view y_min_view y_max_view y_max_view], ...
                     [1 0.9 0.4], 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');
                
                if ~isempty(final_locs)
                     plot(ax, final_locs, final_pks, 'kv', 'MarkerFaceColor', 'y', 'HandleVisibility', 'off', 'Tag', 'PeakMarker');
                    for k = 1:length(final_locs)
                        text(ax, final_locs(k), final_pks(k) + 1.8, sprintf('%.3f Hz', final_locs(k)), ...
                            'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                    end
                end
            end
        end
        
        function [t, d_out] = process_signal(obj, raw)
            d = double(raw(:));
            [b, a] = butter(obj.BUTTER_ORDER, [0.05 1.2] / (obj.FS_TARGET/2), 'bandpass');
            d_out = filtfilt(b, a, d);
            t = (0:length(d_out)-1)' / obj.FS_TARGET;
        end
        
        function cleanAxes(obj, ax)
            cla(ax, 'reset'); 
            grid(ax, 'on');
            delete(findall(ax, 'Type', 'text')); 
            delete(findall(ax, 'Type', 'legend'));
        end
        
        function resetApp(obj)
            obj.FileListBox.Value = {};
            obj.DataStruct = [];
            axesList = [obj.AxTemp1, obj.AxTemp2, obj.AxFreq1, obj.AxFreq2];
            for i = 1:4
                obj.cleanAxes(axesList(i));
            end
            obj.applyTitles();
            xlim(obj.AxFreq1, obj.DISPLAY_FREQ_LIMIT); xlim(obj.AxFreq2, obj.DISPLAY_FREQ_LIMIT);
            obj.loadFileList(); 
        end
    end
end