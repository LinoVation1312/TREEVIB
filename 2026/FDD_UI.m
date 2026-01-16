classdef FDD_UI < handle
    properties (Access = public)
        %% ========== CONFIGURATION (MODIFIABLE ICI) ==========
        MAT_FOLDER = "DATA5HZ";
        FS_TARGET = 5;              
        DEFAULT_FMIN = 0.1;          
        DEFAULT_FMAX = 0.75;         
        DEFAULT_DELTAF = 0.1;        
        DEFAULT_DYN = 10;            
        DEFAULT_NFFT = 1500;         
        BUTTER_ORDER = 2;            
        DISPLAY_FREQ_LIMIT = [0 1];  
        
        % Paramètres de dynamique d'affichage (Y)
        Y_MARGIN_UP = 5;             
        Y_DYNAMIC_RANGE = 30;        
        
        % Couleurs
        ColNord = [0 0.4470 0.7410]; % Bleu
        ColEst  = [0.8500 0.3250 0.0980]; % Orange
        ColSVD  = [0.15 0.15 0.15];  % Anthracite
        
        TreesConfig = struct('Arbre1', {{'Arbre1_Nord', 'Arbre1_Est'}}, ...
                             'Arbre2', {{'Arbre2_Nord', 'Arbre2_Est'}})
    end
    
    properties (Access = private)
        UIFigure, FileListBox
        AxTemp1, AxTemp2
        AxFreq1, AxFreq2
        FMinSpin, FMaxSpin, NFFTEdit
        DeltaFSpin, DynEdit 
        DataStruct
    end

    methods
        function obj = FDD_UI()
            clc; close all; 
            obj.createUI();
            obj.loadFileList();
        end

        function createUI(obj)
            obj.UIFigure = uifigure('Name', 'Analyse Modale FDD - Multi-Arbres', 'Position', [100 100 1300 850]);
            gl = uigridlayout(obj.UIFigure, [1 2]);
            gl.ColumnWidth = {260, '1x'};

            p = uipanel(gl, 'Title', 'Configuration');
            uilabel(p, 'Position', [10 650 200 22], 'Text', 'Mesures disponibles :');
            
            % ListBox modifiée pour gérer Items (texte) et ItemsData (fichier)
            obj.FileListBox = uilistbox(p, 'Position', [10 450 230 180], ...
                'ValueChangedFcn', @(src,event) obj.loadFile());
            
            y = 410;
            uilabel(p, 'Position', [10 y 120 22], 'Text', 'Plage fmin (Hz)');
            obj.FMinSpin = uieditfield(p, 'numeric', 'Position', [140 y 60 22], 'Value', obj.DEFAULT_FMIN);
            y = y-30;
            uilabel(p, 'Position', [10 y 120 22], 'Text', 'Plage fmax (Hz)');
            obj.FMaxSpin = uieditfield(p, 'numeric', 'Position', [140 y 60 22], 'Value', obj.DEFAULT_FMAX);
            y = y-30;
            uilabel(p, 'Position', [10 y 120 22], 'Text', 'Delta f (+/- Hz)');
            obj.DeltaFSpin = uieditfield(p, 'numeric', 'Position', [140 y 60 22], 'Value', obj.DEFAULT_DELTAF);
            y = y-30;
            uilabel(p, 'Position', [10 y 120 22], 'Text', 'Dyn. (dB sous max)');
            obj.DynEdit = uieditfield(p, 'numeric', 'Position', [140 y 60 22], 'Value', obj.DEFAULT_DYN);
            y = y-30;
            uilabel(p, 'Position', [10 y 120 22], 'Text', 'NFFT');
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
            title(obj.AxFreq1, 'Arbre 1 : DSP & SVD'); title(obj.AxFreq2, 'Arbre 2 : DSP & SVD');
            ylabel(obj.AxFreq1, 'dB/Hz'); ylabel(obj.AxFreq2, 'dB/Hz');
            xlabel(obj.AxFreq1, 'Fréquence (Hz)'); xlabel(obj.AxFreq2, 'Fréquence (Hz)');
        end

        %% --- MODIFICATION PRINCIPALE ICI ---
function loadFileList(obj)
            % 1. Vérification du dossier
            if ~exist(obj.MAT_FOLDER, 'dir')
                obj.FileListBox.Items = "Dossier introuvable";
                return; 
            end
            
            % 2. Récupération des fichiers
            files = dir(fullfile(obj.MAT_FOLDER, '*.mat'));
            if isempty(files)
                obj.FileListBox.Items = "Aucun fichier .mat";
                return;
            end
            
            % 3. Tri et initialisation
            [~, idx] = sort({files.name});
            files = files(idx);
            rawNames = string({files.name});
            displayNames = rawNames; 
            
            % 4. Boucle de formatage
            for i = 1:length(rawNames)
                filename = rawNames(i);
                
                % Regex
                tokens = regexp(filename, '(\d{8})_(\d{6})', 'tokens', 'once');
                
                if ~isempty(tokens)
                    try
                        % --- CORRECTION ICI ---
                        % On force la conversion en string pour coller le texte
                        % au lieu d'additionner des matrices
                        dateStr = string(tokens{1}) + string(tokens{2}); 
                        
                        % Conversion DateTime
                        dt = datetime(dateStr, 'InputFormat', 'yyyyMMddHHmmss');
                        
                        % Arrondi minute
                        dt_rounded = dateshift(dt, 'start', 'minute', 'nearest');
                        
                        % Affichage final
                        displayNames(i) = string(dt_rounded, 'dd/MM - HH''h''mm');
                    catch
                        % En cas d'erreur, on garde le nom brut (silencieusement maintenant)
                    end
                end
            end

            % 5. Mise à jour de l'interface
            obj.FileListBox.Items = displayNames;      
            obj.FileListBox.ItemsData = rawNames;      
        end

        function loadFile(obj)
            % Grâce à ItemsData, 'selected' contient le vrai nom de fichier (ex: 2025...mat)
            % même si l'utilisateur a cliqué sur "13/03 - 18h00"
            selected = obj.FileListBox.Value;
            
            if isempty(selected) || strcmp(selected, 'Dossier introuvable'), return; end
            
            try
                fullPath = fullfile(obj.MAT_FOLDER, selected);
                obj.DataStruct = load(fullPath);
                
                % Petit feedback console
                fprintf('Chargement : %s\n', selected);
                obj.plotPreview();
            catch ME
                uialert(obj.UIFigure, ['Erreur chargement : ' ME.message], 'Erreur fichier');
            end
        end
        %% ----------------------------------

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
            
            % Sécurité si les champs n'existent pas dans le MAT
            if ~isfield(obj.DataStruct, names{1}) || ~isfield(obj.DataStruct, names{2})
                text(axT, 0.5, 0.5, 'Données manquantes', 'HorizontalAlignment', 'center');
                return;
            end
            
            [t, dN] = obj.process_signal(obj.DataStruct.(names{1}));
            [~, dE] = obj.process_signal(obj.DataStruct.(names{2}));
            
            % Temporel
            plot(axT, t, dN, 'Color', obj.ColNord, 'DisplayName', 'Nord');
            plot(axT, t, dE, 'Color', obj.ColEst, 'DisplayName', 'Est');
            legend(axT, 'show', 'FontSize', 7);
            
            % DSP (Aperçu)
            nfft = obj.NFFTEdit.Value;
            [Pnn, f] = pwelch(dN, hamming(nfft), nfft/2, nfft, obj.FS_TARGET);
            [Pee, ~] = pwelch(dE, hamming(nfft), nfft/2, nfft, obj.FS_TARGET);
            
            pnn_db = 10*log10(Pnn);
            pee_db = 10*log10(Pee);
            
            plot(axF, f, pnn_db, 'Color', [obj.ColNord 0.4], 'LineWidth', 0.8, 'DisplayName', 'DSP Nord');
            plot(axF, f, pee_db, 'Color', [obj.ColEst 0.4], 'LineWidth', 0.8, 'DisplayName', 'DSP Est');
            
            % --- APPLICATION DE LA DYNAMIQUE DES LE PRE-AFFICHAGE ---
            f_min = obj.FMinSpin.Value; f_max = obj.FMaxSpin.Value;
            in_range = (f >= f_min) & (f <= f_max);
            if any(in_range)
                max_val = max([max(pnn_db(in_range)), max(pee_db(in_range))]);
                y_max_view = max_val + obj.Y_MARGIN_UP;
                y_min_view = y_max_view - obj.Y_DYNAMIC_RANGE;
                ylim(axF, [y_min_view, y_max_view]);
            end
            
            xlim(axF, obj.DISPLAY_FREQ_LIMIT);
            legend(axF, 'show', 'FontSize', 7);
        end

        function runFDD(obj)
            if isempty(obj.DataStruct), return; end
            obj.computeSVD(obj.AxFreq1, 'Arbre1');
            obj.computeSVD(obj.AxFreq2, 'Arbre2');
        end

        function computeSVD(obj, ax, treeKey)
            delete(findall(ax, 'Tag', 'SVD_Plot'));
            delete(findall(ax, 'Type', 'text'));
            delete(findall(ax, 'Type', 'fill'));
            
            names = obj.TreesConfig.(treeKey);
            if ~isfield(obj.DataStruct, names{1}), return; end

            [~, dN] = obj.process_signal(obj.DataStruct.(names{1}));
            [~, dE] = obj.process_signal(obj.DataStruct.(names{2}));
            
            nfft = obj.NFFTEdit.Value;
            [Pnn, f] = pwelch(dN, hamming(nfft), nfft/2, nfft, obj.FS_TARGET);
            [Pee, ~] = pwelch(dE, hamming(nfft), nfft/2, nfft, obj.FS_TARGET);
            [Pne, ~] = cpsd(dN, dE, hamming(nfft), nfft/2, nfft, obj.FS_TARGET);
            
            s1 = zeros(length(f), 1);
            for j = 1:length(f)
                Gyy = [Pnn(j), Pne(j); conj(Pne(j)), Pee(j)];
                [~, S, ~] = svd(Gyy); s1(j) = S(1,1);
            end
            psd_db = 10*log10(s1);
            
            plot(ax, f, psd_db, 'Color', obj.ColSVD, 'LineWidth', 2, 'DisplayName', 'SVD', 'Tag', 'SVD_Plot');

            f_min = obj.FMinSpin.Value; f_max = obj.FMaxSpin.Value;
            df = obj.DeltaFSpin.Value;  dyn_detect = obj.DynEdit.Value;

            in_range = (f >= f_min) & (f <= f_max);
            if any(in_range)
                [max_val, idx_in_sub] = max(psd_db(in_range));
                sub_f = f(in_range);
                f_at_max = sub_f(idx_in_sub);
                search_min = max(f_min, f_at_max - df);
                search_max = min(f_max, f_at_max + df);
                min_h = max_val - dyn_detect;
                
                [pks, locs] = findpeaks(psd_db, f, 'MinPeakHeight', min_h);
                mask = (locs >= search_min) & (locs <= search_max);
                f_res = locs(mask); a_res = pks(mask);

                if ~isempty(f_res)
                    plot(ax, f_res, a_res, 'kv', 'MarkerFaceColor', 'y', 'HandleVisibility', 'off');
                    for k = 1:length(f_res)
                        text(ax, f_res(k), a_res(k) + 1.8, sprintf('%.3f Hz', f_res(k)), ...
                            'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                    end
                end
                
                y_max_view = max_val + obj.Y_MARGIN_UP;
                y_min_view = y_max_view - obj.Y_DYNAMIC_RANGE;
                ylim(ax, [y_min_view, y_max_view]);

                fill(ax, [search_min search_max search_max search_min], [y_min_view y_min_view y_max_view y_max_view], ...
                    [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.1, 'HandleVisibility', 'off');
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
            obj.loadFileList(); % Recharger la liste si jamais le dossier a changé
        end
    end
end