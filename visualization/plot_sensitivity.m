% =========================================================================
% plot_sensitivity.m
% 目的：读取敏感性分析的 .mat 文件，生成高质量柱状图并自动保存
% =========================================================================

out_dir = '../figures/sensitivity';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

exp_dir = '../experiments';
mat_files = dir(fullfile(exp_dir, 'sensitivity_results_para_*.mat'));
if isempty(mat_files)
    error('未找到敏感性分析的 .mat 文件，请先运行 experiments/run_sensitivity_analysis_para.m');
end

[~, idx] = sort([mat_files.datenum], 'descend');
latest_file = fullfile(exp_dir, mat_files(idx(1)).name);
fprintf('正在加载最新数据文件: %s\n', latest_file);
load(latest_file);

Sensitivity_Table = table(config_names', mean_hv, std_hv, mean_util, mean_nrg, mean_size, ...
    'VariableNames', {'Configuration', 'HV_Mean', 'HV_Std', 'Max_Utility', 'Min_Energy', 'ArchiveSize'});
disp('==================== 敏感性分析数据表 ====================');
disp(Sensitivity_Table);

csv_path = fullfile(out_dir, 'Sensitivity_Table.csv');
writetable(Sensitivity_Table, csv_path);
fprintf('表格已保存至: %s\n', csv_path);

fig = figure('Name', '超参数敏感性分析 (HV)', 'Position', [100, 100, 900, 500]);

bar_colors = repmat([0.2 0.6 0.8], length(mean_hv), 1);
bar_colors(1, :) = [0.8 0.2 0.2];

b = bar(mean_hv, 'FaceColor', 'flat');
b.CData = bar_colors;
hold on;
errorbar(mean_hv, std_hv, 'k.', 'LineWidth', 1.5);

set(gca, 'XTick', 1:length(config_names), 'XTickLabel', config_names, 'XTickLabelRotation', 45);
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');
ylabel('Hypervolume (HV) Score', 'FontWeight', 'bold');
title('Sensitivity of MOO Performance to Subpopulation Weights', 'FontWeight', 'bold');
grid on;

ylim_min = min(mean_hv) * 0.95;
ylim_max = max(mean_hv) * 1.05;
if ylim_min < ylim_max
    ylim([ylim_min, ylim_max]);
end

png_path = fullfile(out_dir, 'sensitivity_hv_bar.png');
fig_path = fullfile(out_dir, 'sensitivity_hv_bar.fig');

saveas(fig, png_path);
saveas(fig, fig_path);
fprintf('柱状图已成功保存至: %s\n', out_dir);