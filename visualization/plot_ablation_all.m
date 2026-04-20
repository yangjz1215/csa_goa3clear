function plot_ablation_all(data_file, output_dir)
% plot_ablation_all - 消融实验完整可视化（收敛曲线/HV柱状图/箱型图/统计表）
if nargin < 1 || isempty(data_file)
    data_file = fullfile('..', 'experiments', 'ablation_results_para_Map1_Medium_*.mat');
    files = dir(data_file);
    if ~isempty(files)
        [~, idx] = sort([files.datenum], 'descend');
        data_file = fullfile(files(idx(1)).folder, files(idx(1)).name);
    else
        error('No ablation results file found');
    end
end
if nargin < 2 || isempty(output_dir)
    output_dir = fullfile('..', 'figures', 'ablation');
end
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

load(data_file);

% ✅ 修复点 1：统一变量命名，对齐现在的 4 个变体
variants = {'proposed', 'no_pareto', 'no_subpop', 'no_goa'};
labels = {'Proposed cSA-GOA', 'w/o Pareto Leader', 'w/o Subpop', 'w/o GOA Repulsion'};
n_vars = length(variants); % 动态获取数量，当前为 4

colors = [
    0.8500, 0.3250, 0.0980;
    0.0000, 0.4470, 0.7410;
    0.9290, 0.6940, 0.1250;
    0.4940, 0.1840, 0.5560
];

%% ========== 图1: 收敛曲线 + HV柱状图 ==========
fig1 = figure('Position', [100, 100, 900, 500]);
set(gcf, 'Color', 'w');

subplot(1, 2, 1);
hold on;
for v = 1:n_vars
    variant = variants{v};
    if isfield(results, variant) && isfield(results.(variant), 'convergence_curves')
        curves = results.(variant).convergence_curves;
        if iscell(curves) && ~isempty(curves)
            max_len = max(cellfun(@length, curves));
            avg_curve = zeros(1, max_len);
            n_runs = 0;
            for run = 1:min(length(curves), 30) % 修改为支持30次均值
                if ~isempty(curves{run})
                    c = curves{run};
                    avg_curve(1:length(c)) = avg_curve(1:length(c)) + c;
                    n_runs = n_runs + 1;
                end
            end
            if n_runs > 0
                avg_curve = avg_curve / n_runs;
                lw = 2.0;
                plot(1:length(avg_curve), avg_curve, '-', 'Color', colors(v,:), 'LineWidth', lw, 'DisplayName', labels{v});
            end
        end
    end
end
xlabel('Generations', 'FontWeight', 'bold');
ylabel('Best Fitness', 'FontWeight', 'bold');
title('Convergence Curves (Average of 30 runs)', 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontName', 'Times New Roman', 'FontSize', 9);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
grid on; ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.5; box on;

subplot(1, 2, 2);
hv_means = zeros(1, n_vars);
hv_stds = zeros(1, n_vars);
for v = 1:n_vars
    variant = variants{v};
    if isfield(results, variant) && isfield(results.(variant), 'hv_values')
        vals = results.(variant).hv_values;
        hv_means(v) = nanmean(vals);
        hv_stds(v) = nanstd(vals);
    end
end
% ✅ 修复点 2：将写死的 5 全部替换为动态的 n_vars
bar(1:n_vars, hv_means, 0.6, 'FaceColor', [0.0000, 0.4470, 0.7410], 'EdgeColor', 'k', 'LineWidth', 1.2);
hold on;
errorbar(1:n_vars, hv_means, hv_stds, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'Capsize', 6);
xlabel('Variant', 'FontWeight', 'bold');
ylabel('Hypervolume (HV)', 'FontWeight', 'bold');
title('HV Comparison (Mean ± Std)', 'FontWeight', 'bold');
xticks(1:n_vars);
xticklabels(labels);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
grid on; ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.5; box on;

saveas(fig1, fullfile(output_dir, 'ablation_convergence_hv.fig'));
saveas(fig1, fullfile(output_dir, 'ablation_convergence_hv.png'));
close(fig1);

%% ========== 图2: 双Y轴柱状图 (性能 vs 计算成本) ==========
fig2 = figure('Position', [100, 100, 800, 500]);
set(gcf, 'Color', 'w');

hv_means = zeros(1, n_vars);
iter_means = zeros(1, n_vars);
for i = 1:n_vars
    v = variants{i};
    if isfield(results, v)
        hv_means(i) = results.(v).mean_hv;
        iter_means(i) = mean(results.(v).iter_counts);
    end
end

yyaxis left
b1 = bar(1:n_vars, hv_means, 0.4, 'FaceColor', [0.0000, 0.4470, 0.7410], 'EdgeColor', 'k', 'LineWidth', 1.2);
ylabel('Hypervolume (HV)', 'FontWeight', 'bold', 'Color', [0.0000, 0.4470, 0.7410]);
ylim([min(hv_means)*0.95, max(hv_means)*1.05]);
set(gca, 'ycolor', [0.0000, 0.4470, 0.7410]);

yyaxis right
hold on;
b2 = bar((1:n_vars)+0.4, iter_means, 0.4, 'FaceColor', [0.8500, 0.3250, 0.0980], 'EdgeColor', 'k', 'LineWidth', 1.2);
ylabel('Actual Iterations (Cost)', 'FontWeight', 'bold', 'Color', [0.8500, 0.3250, 0.0980]);
ylim([0, 350]);
set(gca, 'ycolor', [0.8500, 0.3250, 0.0980]);

xticks((1:n_vars) + 0.2);
xticklabels(labels);
legend([b1, b2], {'Performance (HV)', 'Computational Cost (Iter)'}, 'Location', 'northeast');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.2);
grid on; box on;

saveas(fig2, fullfile(output_dir, 'ablation_performance_cost.fig'));
saveas(fig2, fullfile(output_dir, 'ablation_performance_cost.png'));
close(fig2);

%% ========== 图3: 箱型图 (HV分布) ==========
fig3 = figure('Position', [100, 100, 700, 500]);
set(gcf, 'Color', 'w');

all_hv = [];
group_idx = [];
for i = 1:n_vars
    v = variants{i};
    if isfield(results, v) && isfield(results.(v), 'hv_values')
        hv_data = results.(v).hv_values;
        all_hv = [all_hv; hv_data(:)];
        group_idx = [group_idx; i * ones(length(hv_data(:)), 1)];
    end
end

hold on;
for i = 1:n_vars
    idx = (group_idx == i);
    if any(idx)
        b = boxchart(group_idx(idx), all_hv(idx));
        b.BoxFaceColor = colors(i, :);
        b.BoxFaceAlpha = 0.6;
        b.MarkerStyle = 'o';
        b.MarkerColor = [0.2, 0.2, 0.2];
        b.LineWidth = 1.5;
    end
end

xticks(1:n_vars);
xticklabels(labels);
ylabel('Hypervolume (HV)', 'FontWeight', 'bold');
title('HV Distribution (Robustness Analysis - 30 Runs)', 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.2);
grid on; ax = gca; ax.GridLineStyle = '--'; ax.GridAlpha = 0.3; box on;

saveas(fig3, fullfile(output_dir, 'ablation_boxplot.fig'));
saveas(fig3, fullfile(output_dir, 'ablation_boxplot.png'));
close(fig3);

%% ========== 图4: 统计检验表 ==========
fig4 = figure('Position', [100, 100, 700, 300]);
set(gcf, 'Color', 'w');

variant_labels = {'Proposed cSA-GOA', 'w/o Pareto Leader', 'w/o Multi-Subpop', 'w/o GOA Repulsion'};

if ~isfield(results, 'proposed')
    fprintf('Error: proposed variant not found in data\n');
    close(fig4);
    return;
end

proposed_hv = results.proposed.hv_values;

p_values = zeros(n_vars, 1);
h_stats = zeros(n_vars, 1);
mean_improvement = zeros(n_vars, 1);

for v = 1:n_vars
    variant = variants{v};
    if isfield(results, variant) && strcmp(variant, 'proposed')
        p_values(v) = 1.0;
        h_stats(v) = 0;
        mean_improvement(v) = 0;
    elseif isfield(results, variant)
        variant_hv = results.(variant).hv_values;
        [p_values(v), h_stats(v)] = signrank(proposed_hv, variant_hv);
        mean_improvement(v) = (nanmean(proposed_hv) - nanmean(variant_hv)) / nanmean(variant_hv) * 100;
    else
        p_values(v) = NaN; h_stats(v) = NaN; mean_improvement(v) = NaN;
    end
end

sig_labels = cell(n_vars, 1);
for v = 1:n_vars
    if p_values(v) < 0.001
        sig_labels{v} = '***';
    elseif p_values(v) < 0.01
        sig_labels{v} = '**';
    elseif p_values(v) < 0.05
        sig_labels{v} = '*';
    else
        sig_labels{v} = 'ns';
    end
end

col_labels = {'Variant', 'HV (Mean±Std)', 'p-value', 'Sig.', 'Improve (%)'};
table_data = cell(n_vars, 5);
for v = 1:n_vars
    variant = variants{v};
    if isfield(results, variant)
        hv_mean = nanmean(results.(variant).hv_values);
        hv_std = nanstd(results.(variant).hv_values);
        table_data{v, 1} = variant_labels{v};
        table_data{v, 2} = sprintf('%.4f±%.4f', hv_mean, hv_std);
        table_data{v, 3} = sprintf('%.4e', p_values(v));
        table_data{v, 4} = sig_labels{v};
        table_data{v, 5} = sprintf('%.2f%%', mean_improvement(v));
    end
end

uitable('Data', table_data, 'ColumnName', col_labels, ...
    'Units', 'normalized', 'Position', [0.05, 0.1, 0.9, 0.85], ...
    'FontName', 'Times New Roman', 'FontSize', 11);
title('Statistical Tests (Wilcoxon Signed-Rank Test, 30 Runs)', 'FontName', 'Times New Roman', 'FontWeight', 'bold');

saveas(fig4, fullfile(output_dir, 'ablation_statistical_tests.fig'));
saveas(fig4, fullfile(output_dir, 'ablation_statistical_tests.png'));
close(fig4);

fig5 = figure('Position', [100, 100, 800, 600]);
set(gcf, 'Color', 'w');
view(45, 30);
hold on;

marker_styles = {'o', 's', '^', 'd'};
sizes = [60, 50, 50, 50];

for v = 1:n_vars
    variant = variants{v};
    if ~isfield(results, variant) || ~isfield(results.(variant), 'pareto_fronts')
        continue;
    end
    pfs = results.(variant).pareto_fronts;
    if isempty(pfs)
        continue;
    end
    all_u = []; all_l = []; all_e = [];
    for run = 1:length(pfs)
        pf = pfs{run};
        if ~isempty(pf) && size(pf, 2) == 3
            all_u = [all_u; pf(:,1)];
            all_l = [all_l; pf(:,2)];
            all_e = [all_e; pf(:,3)];
        end
    end

    if ~isempty(all_u)
        scatter3(all_l, all_e, all_u, sizes(v), marker_styles{v}, ...
            'filled', 'MarkerFaceColor', colors(v,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1.0, ...
            'DisplayName', labels{v});
    end
end

xlabel('Total Latency (s) \downarrow', 'FontWeight', 'bold');
ylabel('Total Energy (J) \downarrow', 'FontWeight', 'bold');
zlabel('System Utility \uparrow', 'FontWeight', 'bold');
title('3D Pareto Front - Ablation Study', 'FontWeight', 'bold');
legend('Location', 'best', 'FontName', 'Times New Roman', 'FontSize', 9);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
grid on; box on;
view(45, 30);

saveas(fig5, fullfile(output_dir, 'ablation_pareto_3d.fig'));
saveas(fig5, fullfile(output_dir, 'ablation_pareto_3d.png'));
close(fig5);

fig6 = figure('Name', 'Pareto Projections', 'Position', [100, 100, 1500, 450]);
set(gcf, 'Color', 'w');

proj_configs = {
    2, 1, 'Total Latency (s) \downarrow', 'System Utility \uparrow', 1;
    3, 1, 'Total Energy (J) \downarrow', 'System Utility \uparrow', 2;
    2, 3, 'Total Latency (s) \downarrow', 'Total Energy (J) \downarrow', 3
};

for p = 1:3
    subplot(1, 3, proj_configs{p, 5});
    hold on; grid on;

    for v = n_vars:-1:1
        variant = variants{v};
        if ~isfield(results, variant) || ~isfield(results.(variant), 'pareto_fronts'), continue; end

        pfs = results.(variant).pareto_fronts;
        all_data = [];
        for run = 1:length(pfs)
            if ~isempty(pfs{run}), all_data = [all_data; pfs{run}]; end
        end

        if ~isempty(all_data)
            if strcmp(variant, 'proposed')
                scatter(all_data(:, proj_configs{p, 1}), all_data(:, proj_configs{p, 2}), ...
                    40, marker_styles{v}, ...
                    'MarkerEdgeColor', colors(v,:), ...
                    'LineWidth', 1.5, ...
                    'DisplayName', labels{v});
            else
                scatter(all_data(:, proj_configs{p, 1}), all_data(:, proj_configs{p, 2}), ...
                    20, colors(v,:), marker_styles{v}, 'filled', ...
                    'MarkerFaceAlpha', 0.2, ...
                    'MarkerEdgeColor', 'none', ...
                    'DisplayName', labels{v});
            end
        end
    end

    xlabel(proj_configs{p, 3}, 'FontWeight', 'bold');
    ylabel(proj_configs{p, 4}, 'FontWeight', 'bold');
    title(['Projection: ', proj_configs{p, 3}, ' vs ', proj_configs{p, 4}]);

    if p == 1, legend('Location', 'best', 'FontSize', 7); end
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 10, 'LineWidth', 1.1);
    box on;
end

saveas(fig6, fullfile(output_dir, 'ablation_projections_2d.png'));
saveas(fig6, fullfile(output_dir, 'ablation_projections_2d.fig'));
close(fig6);

fprintf('Ablation visualization saved to %s\n', output_dir);
end