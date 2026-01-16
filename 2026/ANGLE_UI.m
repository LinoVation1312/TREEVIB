classdef ANGLE_UI < handle
    % ANGLE_UI_V3 : Interface avec curseurs temporels (Sans Zoom X)
    
    properties (Access = public)
        %% --- CONFIGURATION ---
        MAT_FOLDER = "DATA5HZ"; 
        FS_TARGET = 5;             
        DEFAULT_FS = 5;
        
        % Paramètres de Filtrage
        F_DC = 0.12;                
        FC = 0.6;   
        HP_TAPS = 401;           
        LP_ORDER = 2;
        
        % Couleurs
        ColNord = [0 0.447 0.741];   
        ColEst  = [0.85 0.325 0.098];
        ColSel  = [0.2 0.2 0.2];    
        ColGray = [0.8 0.8 0.8];     
        ColCurs = [0.1 0.1 0.1];    
    end
    
    properties (Access = private)
        % --- Composants UI ---
        UIFigure, MainGrid
        FileListBox
        RangeSlider, EditMin, EditMax
        
        % Labels Résultats
        LblAngle1, LblRatio1
        LblAngle2, LblRatio2
        
        % Axes
        AxTime1, AxOrb1
        AxTime2, AxOrb2
        
        % --- Objets Graphiques ---
        % Arbre 1
        PlotT1_N_All, PlotT1_N_Sel
        PlotT1_E_All, PlotT1_E_Sel
        PlotOrb1_All, PlotOrb1_Sel
        LineMaj1, LineMin1, ZeroLine1
        Curs1_Min, Curs1_Max % NOUVEAU: Curseurs
        
        % Arbre 2
        PlotT2_N_All, PlotT2_N_Sel
        PlotT2_E_All, PlotT2_E_Sel
        PlotOrb2_All, PlotOrb2_Sel
        LineMaj2, LineMin2, ZeroLine2
        Curs2_Min, Curs2_Max % NOUVEAU: Curseurs
        
        % --- Données ---
        CurrentData
        TreesConfig = struct('Arbre1', {{'Arbre1_Nord', 'Arbre1_Est'}}, ...
                             'Arbre2', {{'Arbre2_Nord', 'Arbre2_Est'}})
    end
    
    methods
        %% CONSTRUCTEUR
        function obj = ANGLE_UI()
            close all; clc;
            obj.createUI();
            obj.loadFileList();
        end
        
        %% CREATION DE L'INTERFACE
        function createUI(obj)
            obj.UIFigure = uifigure('Name', 'Analyse Arbres : Mode Curseurs', 'Position', [50 50 1400 900]);
            
            % Layout
            obj.MainGrid = uigridlayout(obj.UIFigure, [1 2]);
            obj.MainGrid.ColumnWidth = {300, '1x'};
            
            % === 1. PANNEAU DE GAUCHE ===
            p = uipanel(obj.MainGrid, 'Title', 'Contrôles');
            p.Layout.Row = 1; p.Layout.Column = 1;
            
            uilabel(p, 'Position', [10 830 200 20], 'Text', 'Fichiers MAT :', 'FontWeight', 'bold');
            obj.FileListBox = uilistbox(p, 'Position', [10 600 280 230], ...
                'ValueChangedFcn', @(src,e) obj.loadFile());
            
            % Sliders
            panelT = uipanel(p, 'Title', 'Sélection Temporelle', 'Position', [10 450 280 140]);
            
            % Slider : Le mouvement déplace les curseurs, le relâchement calcule la PCA
            obj.RangeSlider = uislider(panelT, 'range', 'Position', [20 80 240 3], ...
                'ValueChangedFcn', @(src,e) obj.onSliderReleased(e), ...
                'ValueChangingFcn', @(src,e) obj.onSliderMoving(e)); 
            
            uilabel(panelT, 'Position', [10 40 40 20], 'Text', 'Min:');
            obj.EditMin = uieditfield(panelT, 'numeric', 'Position', [45 40 70 22], ...
                'ValueChangedFcn', @(src,e) obj.onEditChange());
            uilabel(panelT, 'Position', [135 40 40 20], 'Text', 'Max:');
            obj.EditMax = uieditfield(panelT, 'numeric', 'Position', [170 40 70 22], ...
                'ValueChangedFcn', @(src,e) obj.onEditChange());
            uibutton(panelT, 'Text', 'Vue Complète', 'Position', [10 10 260 25], ...
                'ButtonPushedFcn', @(src,e) obj.resetZoom());
            
            % Résultats
            obj.createResultPanel(p, 1, 280);
            obj.createResultPanel(p, 2, 160);
            
            % === 2. ZONE DROITE ===
            gGrid = uigridlayout(obj.MainGrid, [2 2]);
            gGrid.RowHeight = {'1x', '1x'};
            gGrid.ColumnWidth = {'2x', '1x'}; 
            gGrid.RowSpacing = 20; gGrid.ColumnSpacing = 20;
            
            % -- Arbre 1 --
            obj.AxTime1 = uiaxes(gGrid);
            obj.AxTime1.Layout.Row = 1; obj.AxTime1.Layout.Column = 1;
            title(obj.AxTime1, 'Arbre 1 : Vue Globale'); grid(obj.AxTime1, 'on');
            
            obj.AxOrb1 = uiaxes(gGrid);
            obj.AxOrb1.Layout.Row = 1; obj.AxOrb1.Layout.Column = 2;
            title(obj.AxOrb1, 'Trajectoire (Sélection)'); grid(obj.AxOrb1, 'on'); 
            axis(obj.AxOrb1, 'equal'); 
            
            % -- Arbre 2 --
            obj.AxTime2 = uiaxes(gGrid);
            obj.AxTime2.Layout.Row = 2; obj.AxTime2.Layout.Column = 1;
            title(obj.AxTime2, 'Arbre 2 : Vue Globale'); grid(obj.AxTime2, 'on');
            xlabel(obj.AxTime2, 'Temps (s)');
            
            obj.AxOrb2 = uiaxes(gGrid);
            obj.AxOrb2.Layout.Row = 2; obj.AxOrb2.Layout.Column = 2;
            title(obj.AxOrb2, 'Trajectoire (Sélection)'); grid(obj.AxOrb2, 'on'); 
            axis(obj.AxOrb2, 'equal'); 
            
            obj.initPlotObjects();
        end
        
        function createResultPanel(obj, parentPanel, treeID, yPos)
            pRes = uipanel(parentPanel, 'Title', sprintf('Résultats Arbre %d', treeID), ...
                'Position', [10 yPos 280 110]);
            lblA = uilabel(pRes, 'Position', [10 60 260 25], 'Text', 'Direction: -', ...
                'FontWeight', 'bold', 'FontSize', 14);
            lblR = uilabel(pRes, 'Position', [10 30 260 22], 'Text', 'Ratio: -');
            
            if treeID == 1, obj.LblAngle1 = lblA; obj.LblRatio1 = lblR; pRes.BackgroundColor = [0.95 0.98 1];
            else, obj.LblAngle2 = lblA; obj.LblRatio2 = lblR; pRes.BackgroundColor = [1 0.96 0.92]; end
        end
        
        function initPlotObjects(obj)
            % Init Arbre 1
            hold(obj.AxTime1, 'on');
            obj.ZeroLine1 = plot(obj.AxTime1, [NaN NaN], [0 0], 'Color', [0.6 0.6 0.6]);
            % Fond (Gris)
            obj.PlotT1_N_All = plot(obj.AxTime1, NaN, NaN, 'Color', [0.8 0.8 0.8]);
            obj.PlotT1_E_All = plot(obj.AxTime1, NaN, NaN, 'Color', [0.8 0.8 0.8]);
            % Sélection (Couleur)
            obj.PlotT1_N_Sel = plot(obj.AxTime1, NaN, NaN, 'Color', obj.ColNord, 'LineWidth', 1.5, 'DisplayName', 'Nord');
            obj.PlotT1_E_Sel = plot(obj.AxTime1, NaN, NaN, 'Color', obj.ColEst, 'LineWidth', 1.5, 'DisplayName', 'Est');
            % Curseurs (xline)
            obj.Curs1_Min = xline(obj.AxTime1, NaN, '--k', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
            obj.Curs1_Max = xline(obj.AxTime1, NaN, '--k', 'LineWidth', 1.5);
            legend(obj.AxTime1, [obj.PlotT1_N_Sel, obj.PlotT1_E_Sel], 'Location', 'northwest', 'AutoUpdate', 'off');
            
            hold(obj.AxOrb1, 'on');
            obj.PlotOrb1_All = plot(obj.AxOrb1, NaN, NaN, '.', 'Color', obj.ColGray, 'MarkerSize', 4);
            obj.LineMaj1 = plot(obj.AxOrb1, NaN, NaN, 'r-', 'LineWidth', 2);
            obj.LineMin1 = plot(obj.AxOrb1, NaN, NaN, 'g-', 'LineWidth', 2);
            obj.PlotOrb1_Sel = plot(obj.AxOrb1, NaN, NaN, '-', 'Color', obj.ColSel, 'LineWidth', 1.0);

            % Init Arbre 2
            hold(obj.AxTime2, 'on');
            obj.ZeroLine2 = plot(obj.AxTime2, [NaN NaN], [0 0], 'Color', [0.6 0.6 0.6]);
            obj.PlotT2_N_All = plot(obj.AxTime2, NaN, NaN, 'Color', [0.8 0.8 0.8]);
            obj.PlotT2_E_All = plot(obj.AxTime2, NaN, NaN, 'Color', [0.8 0.8 0.8]);
            obj.PlotT2_N_Sel = plot(obj.AxTime2, NaN, NaN, 'Color', obj.ColNord, 'LineWidth', 1.5);
            obj.PlotT2_E_Sel = plot(obj.AxTime2, NaN, NaN, 'Color', obj.ColEst, 'LineWidth', 1.5);
            % Curseurs
            obj.Curs2_Min = xline(obj.AxTime2, NaN, '--k', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
            obj.Curs2_Max = xline(obj.AxTime2, NaN, '--k', 'LineWidth', 1.5);
            
            hold(obj.AxOrb2, 'on');
            obj.PlotOrb2_All = plot(obj.AxOrb2, NaN, NaN, '.', 'Color', obj.ColGray, 'MarkerSize', 4);
            obj.LineMaj2 = plot(obj.AxOrb2, NaN, NaN, 'r-', 'LineWidth', 2);
            obj.LineMin2 = plot(obj.AxOrb2, NaN, NaN, 'g-', 'LineWidth', 2);
            obj.PlotOrb2_Sel = plot(obj.AxOrb2, NaN, NaN, '-', 'Color', obj.ColSel, 'LineWidth', 1.0);
        end
        
        %% CHARGEMENT
        function loadFileList(obj)
            % 1. Vérification du dossier
            if ~exist(obj.MAT_FOLDER, 'dir')
                obj.FileListBox.Items = "Simulation (Dossier introuvable)";
                obj.FileListBox.ItemsData = "Simulation";
                return; 
            end
            
            % 2. Récupération des fichiers
            files = dir(fullfile(obj.MAT_FOLDER, '*.mat'));
            if isempty(files)
                obj.FileListBox.Items = "Simulation (Aucun fichier)";
                obj.FileListBox.ItemsData = "Simulation";
                return;
            end
            
            % 3. Tri et préparation
            [~, idx] = sort({files.name});
            files = files(idx);
            
            rawNames = string({files.name});
            displayNames = rawNames; 
            
            % 4. Boucle de formatage robuste (Correction Dimension Array)
            for i = 1:length(rawNames)
                filename = rawNames(i);
                
                % Regex : cherche 8 chiffres + underscore + 6 chiffres
                tokens = regexp(filename, '(\d{8})_(\d{6})', 'tokens', 'once');
                
                if ~isempty(tokens)
                    try
                        % Concaténation forcée en string pour éviter l'erreur de dimension
                        dateStr = string(tokens{1}) + string(tokens{2});
                        
                        % Conversion DateTime
                        dt = datetime(dateStr, 'InputFormat', 'yyyyMMddHHmmss');
                        
                        % Arrondi à la minute
                        dt_rounded = dateshift(dt, 'start', 'minute', 'nearest');
                        
                        % Affichage lisible
                        displayNames(i) = string(dt_rounded, 'dd/MM - HH''h''mm');
                    catch
                        % En cas d'erreur, on garde le nom brut
                    end
                end
            end
            
            % 5. Mise à jour Interface (Items = Visuel, ItemsData = Vrai nom)
            obj.FileListBox.Items = displayNames;
            obj.FileListBox.ItemsData = rawNames;
        end

        function loadFile(obj)
            filename = obj.FileListBox.Value;
            
            if isempty(filename), return; end
            
            % Mise à jour visuelle rapide
            obj.LblAngle1.Text = '...'; obj.LblAngle2.Text = '...'; drawnow;
            
            % Gestion du cas Simulation vs Fichier Réel
            if strcmp(filename, "Simulation")
                obj.generateSimulation();
            else
                fullpath = fullfile(obj.MAT_FOLDER, filename);
                try
                    raw = load(fullpath);
                    obj.CurrentData.RawStruct = raw;
                    
                    % Détection FS (Fréquence échantillonnage)
                    fs = obj.DEFAULT_FS;
                    if isfield(raw, 'fs'), fs = double(raw.fs); 
                    elseif isfield(raw, 'FS'), fs = double(raw.FS); end
                    
                    obj.CurrentData.FsOrigin = fs;
                    obj.CurrentData.IsSim = false;
                    obj.processAllData();
                    
                    % Petit log console pour confirmer le chargement
                    fprintf('Chargé : %s\n', filename);
                    
                catch ME
                    uialert(obj.UIFigure, ME.message, 'Erreur de lecture');
                end
            end
        end
        
        function processAllData(obj)
            fs0 = obj.CurrentData.FsOrigin;
            if obj.CurrentData.IsSim
                t = obj.CurrentData.t;
                d1N = obj.CurrentData.A1N; d1E = obj.CurrentData.A1E;
                d2N = obj.CurrentData.A2N; d2E = obj.CurrentData.A2E;
            else
                S = obj.CurrentData.RawStruct;
                try, [t, d1N] = obj.process_signal(S.(obj.TreesConfig.Arbre1{1}), fs0);
                     [~, d1E] = obj.process_signal(S.(obj.TreesConfig.Arbre1{2}), fs0);
                catch, d1N=[]; d1E=[]; end
                try, [~, d2N] = obj.process_signal(S.(obj.TreesConfig.Arbre2{1}), fs0);
                     [~, d2E] = obj.process_signal(S.(obj.TreesConfig.Arbre2{2}), fs0);
                catch, d2N=[]; d2E=[]; end
            end
            
            obj.CurrentData.t_proc = t;
            obj.CurrentData.A1N_proc = d1N; obj.CurrentData.A1E_proc = d1E;
            obj.CurrentData.A2N_proc = d2N; obj.CurrentData.A2E_proc = d2E;
            
            if isempty(t), return; end
            
            % Reset Slider
            obj.RangeSlider.Limits = [t(1), t(end)];
            obj.RangeSlider.Value = [t(1), t(end)];
            obj.EditMin.Value = t(1); obj.EditMax.Value = t(end);
            
            % Fixer les limites des axes sur TOUT le fichier (Pas de zoom)
            xlim(obj.AxTime1, [t(1) t(end)]);
            xlim(obj.AxTime2, [t(1) t(end)]);
            
            % Affichage Données de Fond (Gris)
            obj.updateBackground(t, d1N, d1E, obj.PlotT1_N_All, obj.PlotT1_E_All, obj.PlotOrb1_All, obj.ZeroLine1);
            obj.updateBackground(t, d2N, d2E, obj.PlotT2_N_All, obj.PlotT2_E_All, obj.PlotOrb2_All, obj.ZeroLine2);
            
            obj.updateVizFull([t(1), t(end)]);
        end
        
        function updateBackground(obj, t, N, E, hN, hE, hOrb, hZero)
            if isempty(N), return; end
            set(hN, 'XData', t, 'YData', N);
            set(hE, 'XData', t, 'YData', E);
            set(hOrb, 'XData', E, 'YData', N);
            set(hZero, 'XData', [t(1) t(end)]);
            
            % Auto Zoom Y Global pour voir tout le signal verticalement
            mx = max(max(abs(N)), max(abs(E))) * 1.1;
            if mx==0, mx=0.1; end
            if hN == obj.PlotT1_N_All
                ylim(obj.AxTime1, [-mx mx]);
            else
                ylim(obj.AxTime2, [-mx mx]);
            end
        end
        
        %% GESTION SLIDER & CURSEURS
        function onSliderMoving(obj, event)
            % Appelé pendant le glissement : Met à jour les curseurs et champs texte uniquement
            if isprop(event, 'Value'), val = event.Value;
            elseif isprop(event, 'NewValue'), val = event.NewValue; 
            else, return; end
            
            obj.EditMin.Value = val(1); 
            obj.EditMax.Value = val(2);
            
            % Déplace visuellement les lignes sans recalculer la PCA
            obj.Curs1_Min.Value = val(1); obj.Curs1_Max.Value = val(2);
            obj.Curs2_Min.Value = val(1); obj.Curs2_Max.Value = val(2);
        end
        
        function onSliderReleased(obj, event)
            % Appelé au relâchement : Calcul PCA et mise à jour des tracés colorés
            obj.updateVizFull(event.Value);
        end
        
        function onEditChange(obj)
            vMin = obj.EditMin.Value; vMax = obj.EditMax.Value;
            lims = obj.RangeSlider.Limits;
            vMin = max(lims(1), min(vMin, vMax)); 
            vMax = min(lims(2), max(vMin, vMax));
            obj.RangeSlider.Value = [vMin, vMax];
            obj.updateVizFull([vMin, vMax]);
        end
        
        function resetZoom(obj)
             if isempty(obj.CurrentData), return; end
             lims = obj.RangeSlider.Limits;
             obj.RangeSlider.Value = lims;
             obj.updateVizFull(lims);
        end
        
        %% COEUR GRAPHIQUE
        function updateVizFull(obj, rangeVals)
            if isempty(obj.CurrentData) || ~isfield(obj.CurrentData, 't_proc'), return; end
            
            t = obj.CurrentData.t_proc;
            tmin = rangeVals(1); tmax = rangeVals(2);
            obj.EditMin.Value = tmin; obj.EditMax.Value = tmax;
            
            % Mise à jour position Curseurs
            obj.Curs1_Min.Value = tmin; obj.Curs1_Max.Value = tmax;
            obj.Curs2_Min.Value = tmin; obj.Curs2_Max.Value = tmax;
            
            mask = (t >= tmin) & (t <= tmax);
            
            % --- ARBRE 1 ---
            obj.updateTreeViz(t(mask), mask, ...
                obj.CurrentData.A1N_proc, obj.CurrentData.A1E_proc, ...
                obj.PlotT1_N_Sel, obj.PlotT1_E_Sel, obj.PlotOrb1_Sel, ...
                obj.LineMaj1, obj.LineMin1, obj.LblAngle1, obj.LblRatio1);
                
            % --- ARBRE 2 ---
            obj.updateTreeViz(t(mask), mask, ...
                obj.CurrentData.A2N_proc, obj.CurrentData.A2E_proc, ...
                obj.PlotT2_N_Sel, obj.PlotT2_E_Sel, obj.PlotOrb2_Sel, ...
                obj.LineMaj2, obj.LineMin2, obj.LblAngle2, obj.LblRatio2);
        end
        
        function updateTreeViz(obj, t_sub, full_mask, N_full, E_full, ...
                               hT_N, hT_E, hOrb, hLMaj, hLMin, lblAng, lblRat)
            if isempty(N_full) || sum(full_mask) < 2
                set(hT_N, 'XData', NaN, 'YData', NaN); 
                set(hT_E, 'XData', NaN, 'YData', NaN);
                return; 
            end
            
            % Récupération données sélectionnées
            N = N_full(full_mask); 
            E = E_full(full_mask);
            
            % Mettre à jour la surbrillance (le tracé coloré par dessus le gris)
            set(hT_N, 'XData', t_sub, 'YData', N);
            set(hT_E, 'XData', t_sub, 'YData', E);
            
            % Update Orbite (Zoom sur la partie sélectionnée seulement)
            set(hOrb, 'XData', E, 'YData', N);
            
            % PCA
            [ang, rat, LMaj, LMin] = obj.computePCA(E, N);
            lblAng.Text = sprintf('Direction: %.1f°', ang);
            lblRat.Text = sprintf('Ratio: %.2f', rat);
            
            set(hLMaj, 'XData', LMaj(:,1), 'YData', LMaj(:,2));
            set(hLMin, 'XData', LMin(:,1), 'YData', LMin(:,2));
        end
        
        function [angle, ratio, LineMaj, LineMin] = computePCA(~, E, N)
            if length(E) < 5
                angle=NaN; ratio=NaN; LineMaj=[NaN NaN]; LineMin=[NaN NaN]; return;
            end
            mE = mean(E); mN = mean(N);
            C = cov(E - mE, N - mN);
            [V, D] = eig(C);
            [eig_vals, sort_idx] = sort(diag(D), 'descend');
            eig_vecs = V(:, sort_idx);
            
            scale = 2.5; 
            s_maj = scale * sqrt(eig_vals(1)); s_min = scale * sqrt(eig_vals(2));
            v_maj = eig_vecs(:,1); v_min = eig_vecs(:,2);
            
            LineMaj = [mE - v_maj(1)*s_maj, mN - v_maj(2)*s_maj; mE + v_maj(1)*s_maj, mN + v_maj(2)*s_maj];
            LineMin = [mE - v_min(1)*s_min, mN - v_min(2)*s_min; mE + v_min(1)*s_min, mN + v_min(2)*s_min];
                   
            angle = atan2d(v_maj(2), v_maj(1));
            if eig_vals(1) > 0, ratio = sqrt(eig_vals(2)/eig_vals(1)); else, ratio = 0; end
        end
        
        %% PROCESS & SIMU
        function [t, d_out] = process_signal(obj, raw_data, fs_in)
            raw_data = double(raw_data(:)) - mean(double(raw_data(:)));
            N_orig = length(raw_data);
            if fs_in ~= obj.FS_TARGET, d_res = resample(raw_data, round(N_orig*obj.FS_TARGET/fs_in), N_orig);
            else, d_res = raw_data; end
            
            nyq = obj.FS_TARGET/2;
            if obj.F_DC > 0
                taps = obj.HP_TAPS;
                if length(d_res) < taps*1.5, taps = floor(length(d_res)/3); end
                if mod(taps,2)==0, taps=taps-1; end
                if taps > 6, b_hp = fir1(taps-1, obj.F_DC/nyq, 'high'); d_hp = filtfilt(b_hp, 1, d_res);
                else, d_hp = d_res; end
            else, d_hp = d_res; end
            
            [b_lp, a_lp] = butter(obj.LP_ORDER, obj.FC/nyq, 'low');
            d_out = filtfilt(b_lp, a_lp, d_hp);
            t = (0:length(d_out)-1)' / obj.FS_TARGET;
        end

    end
end