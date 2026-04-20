%% CEC 2017 Results Processing and Visualization Script (SCI Paper Standard)
% Function: Generate CSV tables, filled boxplots, and log convergence curves
clear; clc; close all;

filename = 'results_cec2017_D30.mat';
if ~exist(filename, 'file')
    error('Cannot find data file: %s', filename);
end
load(filename);

original_dir = pwd;
cec_core_path = fullfile(original_dir, 'cec2017_data', 'cec2017');
if exist(cec_core_path, 'dir')
    addpath(cec_core_path);
else
    error('CEC data path not found: %s', cec_core_path);
end

func_indices = find(~cellfun(@isempty, {results.mean}));

fprintf('Generating CSV summary table...\n');
csv_file = 'CEC2017_Statistical_Results.csv';
fid = fopen(csv_file, 'w');
fprintf(fid, 'Function,Mean,Std,Best,Worst,Median\n');
for i = 1:length(func_indices)
    f = func_indices(i);
    raw_data = results(f).raw;
    fprintf(fid, 'F%d,%.4e,%.4e,%.4e,%.4e,%.4e\n', ...
        f, results(f).mean, results(f).std, min(raw_data), max(raw_data), median(raw_data));
end
fclose(fid);

fprintf('Generating SCI-standard boxplots...\n');
plot_funcs = [1, 5, 10, 15, 20, 25];
figure('Color', 'w', 'Units', 'pixels', 'Position', [100, 100, 1000, 600]);

faceColor = [0.2, 0.4, 0.7];

for i = 1:length(plot_funcs)
    f = plot_funcs(i);
    subplot(2, 3, i);

    h = boxplot(results(f).raw, 'Labels', {['F', num2str(f)]}, ...
        'Symbol', 'r+', 'OutlierSize', 4, 'Widths', 0.5);

    obj = findobj(gca, 'Tag', 'Box');
    for j = 1:length(obj)
        patch(get(obj(j), 'XData'), get(obj(j), 'YData'), faceColor, 'FaceAlpha', 0.5);
    end

    set(findobj(gca, 'Tag', 'Median'), 'Color', 'r', 'LineWidth', 2);

    ylabel('Error Value');
    grid on;
    set(gca, 'FontSize', 10, 'FontName', 'Arial');
    title(['Function F', num2str(f)], 'FontSize', 12);
end
saveas(gcf, 'CEC2017_Boxplots_Color.png');

fprintf('Generating convergence curves...\n');
figure('Color', 'w', 'Units', 'pixels', 'Position', [150, 150, 1000, 600]);
colors = lines(length(plot_funcs));

for i = 1:length(plot_funcs)
    f = plot_funcs(i);
    fprintf('  Capturing curve for F%d...\n', f);

    addpath(original_dir);
    cd(cec_core_path);
    Dim = 30; max_fes = 300000; Lb = -100 * ones(1, Dim); Ub = 100 * ones(1, Dim);
    [~, ~, curve] = main_cec(f, Dim, max_fes, Lb, Ub);
    cd(original_dir);

    subplot(2, 3, i);
    semilogy(curve, 'LineWidth', 2.5, 'Color', colors(i,:));

    xlabel('FES');
    ylabel('Best Fitness');
    title(['F', num2str(f), ' Convergence'], 'FontSize', 12);
    grid on;
    set(gca, 'FontSize', 10, 'FontName', 'Arial');
end
saveas(gcf, 'CEC2017_Convergence_Curves_New.png');

disp('All paper materials generated successfully!');