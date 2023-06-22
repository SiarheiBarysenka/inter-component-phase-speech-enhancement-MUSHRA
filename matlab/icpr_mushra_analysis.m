function icpr_mushra_analysis
%   ICPR_MUSHRA_ANALYSIS conducts statistical analysis of MUSHRA data
%   using the non-parametric Friedman test. The data is collected from
%   14 listeners exclusively for the work reported in [1].

%   This MATLAB script is utilized in [1] to generate Figure 9,
%   Tables III, and IV.

%   References:
%      [1] S.Y. Barysenka and V.I. Vorobiov, “SNR-based inter-component
%          phase estimation using bi-phase prior statistics for single-
%          channel speech enhancement,” IEEE/ACM Transactions on Audio,
%          Speech and Language Processing, vol. 31, pp. 2365-2381, 2023.
%          DOI: https://doi.org/10.1109/TASLP.2023.3284514

%   Copyright 2023 Siarhei Y. Barysenka

    clear variables; clc; close all;

    set(0, 'DefaultAxesFontSize', 14);
    set(0, 'DefaultTextFontSize', 14);
    set(0, 'DefaultTextInterpreter', 'latex');
    set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
    set(0, 'DefaultLegendInterpreter', 'latex');
    
    %% Load MUSHRA data
    %
    % Loads the data into the following variables:
    %
    %   - conditions (1x7 string):
    %         Names of all systems evaluated in the test:
    %         ["Noisy", "SE+BVM", "BH+BLW", "MMSE-LSA", "MMSE-LSA+SE+BVM", "MMSE-LSA+BH+BLW", "Clean"].
    %
    %   - stimuli (1x6 string):
    %        Descriptions of the noise environments evaluated in the test,
    %        in the format of "<noise_title>-<SNR>":
    %        ["Pink-5", "Pink-10", "Factory-5", "Factory-10", "Babble-5", "Babble-10"].
    %
    %   - data (7x6x14 double):
    %        MUSHRA scores for 14 listeners for each of the 7 systems and 6 noise environments.
    
    load("icpr_mushra_data.mat");
    
    % Convert 3D-array [7 x 6 x 14] into 2D array [7 x (6 * 14)]
    data_for_stats = reshape(data, size(data, 1), size(data, 2) * size(data, 3));
    
    %% Present the table of the mean MUSHRA scores
    mean_values_table = array2table(conditions', "VariableNames", "Condition");
    mean_values_table.("Mean") = mean(data_for_stats, 2);
    disp(mean_values_table);
    
    %% Calculate Friedman statistics
    [~, ~, stats] = friedman(data_for_stats', 1);
    
    %% Show multiple comparisons
    figure();
    alpha = 0.05;
    confidence_interval = (1 - alpha) * 100;
    
    % If 'CType' ('CriticalValueType') parameter is not provided in "multcompare",
    % MATLAB defaults to Tukey-Kramer correction of p-values.
    [comparison, ~, ~, ~] = multcompare(stats, 'Alpha', alpha);
    
    % Decoration of the figure appearance and axes
    figure_handle = gcf;
    set(figure_handle, 'units', 'points', 'position', [0, 0, 460, 180]);
    figure_handle.Children.YTickLabel = flip(conditions);
    figure_handle.Children.XTick = 1 : numel(conditions);
    figure_handle.Children.XLim = [1 (numel(conditions) + 1)];
    figure_handle.Children.Title.String = "";
    figure_handle.Children.XLabel.String = ['Mean Rank' ' (' num2str(confidence_interval, '%.0f') '\% Confidence)'];
    
    % Decoration of lines and markers
    number_of_lines = numel(figure_handle.Children.Children);
    for i = 1 : number_of_lines
        line = figure_handle.Children.Children(i);
        if line.LineStyle == '-'
            % Adjust style of lines
            line.LineWidth = 2;
            line.Marker = 'none';
            line.MarkerSize = 2;
            line.Color = 'black';
        end
        if line.Marker == 'o'
            % Adjust style of markers
            line.Marker = '.';
            line.MarkerSize = 20;
            line.Color = 'black';
        end
    end
    
    % Decoration of grid
    set(gca, 'YGrid', 'on', 'XGrid', 'on');
    axes_handle = gca;
    axes_handle.LineWidth = 1;
    axes_handle.GridLineStyle = '--';
    axes_handle.GridColor = 'black';
    axes_handle.GridAlpha = 0.25;
    
    %% Present table with p-values and effect sizes
    condition_A_label = "Condition A";
    condition_B_label = "Condition B";
    diff_table = array2table(comparison, "VariableNames", [condition_A_label, condition_B_label, "Lower Limit", "A-B", "Upper Limit", "P-value"]);
    
    group_A_indices = diff_table.(condition_A_label);
    group_B_indices = diff_table.(condition_B_label);
    
    diff_table.(condition_A_label) = arrayfun(@(x) conditions(x), group_A_indices);
    diff_table.(condition_B_label) = arrayfun(@(x) conditions(x), group_B_indices);
    
    % Calculate Cliff's Delta
    cliffs_deltas = arrayfun( ...
        @(index) ...
            cliffs_delta( ...
                data_for_stats(group_A_indices(index), :), ...
                data_for_stats(group_B_indices(index), :) ...
            ), ...
        1 : numel(group_A_indices) ...
    );

    % Map Cliff's Deltas to effect sizes
    effect_sizes = arrayfun(@(x) effect_size(x), cliffs_deltas);
    
    diff_table.("Abs of Cliff's Delta") = abs(cliffs_deltas');
    diff_table.("Effect Size") = effect_sizes';
    disp(diff_table);
end

function result = cliffs_delta(x, y)
    m = length(x);
    n = length(y);

    count = 0;
    for i = 1 : m
        for j = 1 : n
            if (x(i) > y(j))
                count = count + 1;
            elseif (x(i) < y(j))
                count = count - 1;
            end
        end
    end

    result = count / (m * n);
end

function result = effect_size(cliffs_delta)
    delta = abs(cliffs_delta);
    
    if (delta < 0.11)
        result = "";
    elseif (0.11 <= delta && delta < 0.28)
        result = "small";
    elseif (0.28 <= delta && delta < 0.43)
        result = "medium";
    else
        result = "large";
    end
end