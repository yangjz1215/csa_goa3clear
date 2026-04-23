function plot_comparison_all(data_file, output_dir)
% plot_comparison_all - 对比实验完整可视化（收敛曲线/Pareto散点/箱型图）
% 输出至 figures/comparison 文件夹
if nargin < 1 || isempty(data_file)
    data_file = fullfile('..', 'experiments', 'comparison_results_para_Map1_Medium_*.mat');
    files = dir(data_file);
    if ~isempty(files)
        [~, idx] = sort([files.datenum], 'descend');
        data_file = fullfile(files(idx(1)).folder, files(idx(1)).name);
    else
        error('No comparison results file found');
    end
end
if nargin < 2 || isempty(output_dir)
    output_dir = fullfile('..', 'figures', 'comparison');
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

load(data_file);

fig1 = figure('Position', [100, 100, 900, 500]);
set(gcf, 'Color', 'w');

colors = [
    0.8500, 0.3250, 0.0980;
    0.0000, 0.4470, 0.7410;
    0.9290, 0.6940, 0.1250;
    0.4940, 0.1840, 0.5560;
    0.4660, 0.6740, 0.1880;
    0.6350, 0.0780, 0.1840
];

algorithms = {'cSA_GOA', 'PSO', 'GA', 'GOA', 'cSA', 'GWO'};
labels = {'cSA-GOA (Proposed)', 'PSO', 'GA', 'GOA', 'cSA', 'GWO'};
marker_styles = {'*', 'o', 's', '^', 'd', 'v'};
sizes = [80, 40, 40, 40, 40, 40];
n_algs = length(algorithms);

subplot(1, 2, 1);
hold on;
for a = 1:length(algorithms)
    alg = algorithms{a};
    if isfield(results, alg) && isfield(results.(alg), 'convergence_curves')
        curves = results.(alg).convergence_curves;
        if iscell(curves) && ~isempty(curves)
            max_len = max(cellfun(@length, curves));
            avg_curve = zeros(1, max_len);
            n_runs = 0;
            for run = 1:min(length(curves), 30)
                if ~isempty(curves{run})
                    c = curves{run};
                    c_min = min(c);
                    c_max = max(c);
                    if c_max > c_min
                        c_norm = (c - c_min) / (c_max - c_min);
                    else
                        c_norm = zeros(size(c));
                    end
                    avg_curve(1:length(c_norm)) = avg_curve(1:length(c_norm)) + c_norm;
                    n_runs = n_runs + 1;
                end
            end
            if n_runs > 0
                avg_curve = avg_curve / n_runs;
                if a == 1
                    lw = 2.5;
                else
                    lw = 1.8;
                end
                plot(1:length(avg_curve), avg_curve, '-', 'Color', colors(a,:), 'LineWidth', lw, 'DisplayName', labels{a});
            end
        end
    end
end
xlabel('Generations', 'FontWeight', 'bold');
ylabel('Normalized Fitness (0-1) \uparrow', 'FontWeight', 'bold');
title('Normalized Convergence', 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontName', 'Times New Roman', 'FontSize', 9);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
ylim([0, 1.05]);
grid on;
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.5;
box on;

subplot(1, 2, 2);
hold on;
hv_means = zeros(1, length(algorithms));
hv_stds = zeros(1, length(algorithms));
for a = 1:length(algorithms)
    alg = algorithms{a};
    if isfield(results, alg) && isfield(results.(alg), 'hv_values')
        vals = results.(alg).hv_values;
        hv_means(a) = mean(vals);
        hv_stds(a) = std(vals);
    end
end
errorbar(1:length(algorithms), hv_means, hv_stds, '-o', 'Color', [0.0000, 0.4470, 0.7410], ...
    'LineWidth', 2.0, 'MarkerSize', 8, 'MarkerFaceColor', [0.0000, 0.4470, 0.7410], 'Capsize', 6);
xlabel('Algorithm', 'FontWeight', 'bold');
ylabel('Hypervolume (HV)', 'FontWeight', 'bold');
title('HV Comparison (Mean ± Std)', 'FontWeight', 'bold');
xticks(1:length(algorithms));
xticklabels(labels);
xtickangle(45);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10, 'LineWidth', 1.2);
grid on;
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.5;
box on;

saveas(fig1, fullfile(output_dir, 'comparison_convergence.fig'));
saveas(fig1, fullfile(output_dir, 'comparison_convergence.png'));
close(fig1);

fig2 = figure('Position', [100, 100, 800, 600]);
set(gcf, 'Color', 'w');
hold on;
grid on;
view(45, 30);

sz = 60;
for a = 1:length(algorithms)
    alg = algorithms{a};
    if isfield(results, alg) && isfield(results.(alg), 'pareto_fronts')
        pfs = results.(alg).pareto_fronts;

        best_run = 1; max_size = 0;
        for run = 1:length(pfs)
            if ~isempty(pfs{run}) && size(pfs{run},1) > max_size
                max_size = size(pfs{run},1);
                best_run = run;
            end
        end

        best_pf = pfs{best_run};
        if ~isempty(best_pf) && size(best_pf,2) == 3
            scatter3(best_pf(:,1), best_pf(:,2), best_pf(:,3), sz, ...
                'MarkerEdgeColor', colors(a,:), 'MarkerFaceColor', colors(a,:), ...
                'MarkerFaceAlpha', 0.6, 'DisplayName', labels{a});
        end
    end
end

xlabel('System Utility (Priority Sum)', 'FontWeight', 'bold');
ylabel('Total Latency (s)', 'FontWeight', 'bold');
zlabel('Total Energy (J)', 'FontWeight', 'bold');
title('3D Pareto Front Comparison', 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontName', 'Times New Roman', 'FontSize', 10);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
box on;

saveas(fig2, fullfile(output_dir, 'comparison_pareto_3D.fig'));
saveas(fig2, fullfile(output_dir, 'comparison_pareto_3D.png'));
close(fig2);

fig3 = figure('Position', [100, 100, 700, 500]);
set(gcf, 'Color', 'w');

hv_all = [];
group_idx = [];
for a = 1:length(algorithms)
    alg = algorithms{a};
    if isfield(results, alg) && isfield(results.(alg), 'hv_values')
        vals = results.(alg).hv_values;
        hv_all = [hv_all; vals(:)];
        group_idx = [group_idx; a * ones(length(vals(:)), 1)];
    end
end

box_colors = [
    0.8500, 0.3250, 0.0980;
    0.0000, 0.4470, 0.7410;
    0.9290, 0.6940, 0.1250;
    0.4940, 0.1840, 0.5560;
    0.4660, 0.6740, 0.1880;
    0.6350, 0.0780, 0.1840
];

hold on;
for a = 1:length(algorithms)
    idx = (group_idx == a);
    if any(idx)
        b = boxchart(group_idx(idx), hv_all(idx));
        b.BoxFaceColor = box_colors(a, :);
        b.BoxFaceAlpha = 0.6;
        b.MarkerStyle = 'o';
        b.MarkerColor = [0.2, 0.2, 0.2];
        b.LineWidth = 1.5;
    end
end

xticks(1:length(algorithms));
xticklabels(labels);
ylabel('Hypervolume (HV)', 'FontWeight', 'bold');
title('HV Distribution (30 Independent Runs)', 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
grid on;
ax = gca; ax.GridLineStyle = '--'; ax.GridAlpha = 0.3;
box on;

saveas(fig3, fullfile(output_dir, 'comparison_boxplot.fig'));
saveas(fig3, fullfile(output_dir, 'comparison_boxplot.png'));
close(fig3);

fprintf('Comparison charts saved to %s\n', output_dir);
fprintf('  - comparison_convergence.fig/png (收敛曲线 + HV误差棒)\n');
fprintf('  - comparison_pareto_3D.fig/png (3D Pareto 前沿)\n');
fprintf('  - comparison_boxplot.fig/png (HV分布箱型图)\n');

fig4 = figure('Name', 'Pareto Projections', 'Position', [100, 100, 1500, 450]);
set(gcf, 'Color', 'w');

proj_configs = {
    2, 1, 'Total Latency (s) \downarrow', 'System Utility \uparrow', 1;
    3, 1, 'Total Energy (J) \downarrow', 'System Utility \uparrow', 2;
    2, 3, 'Total Latency (s) \downarrow', 'Total Energy (J) \downarrow', 3
};

for p = 1:3
    subplot(1, 3, proj_configs{p, 5});
    hold on; grid on;

    for a = 1:length(algorithms)
        alg = algorithms{a};
        if ~isfield(results, alg) || ~isfield(results.(alg), 'pareto_fronts'), continue; end

        pfs = results.(alg).pareto_fronts;
        all_data = [];
        for run = 1:length(pfs)
            if ~isempty(pfs{run}), all_data = [all_data; pfs{run}]; end
        end

        if ~isempty(all_data)
            if strcmp(alg, 'cSA_GOA')
                scatter(all_data(:, proj_configs{p, 1}), all_data(:, proj_configs{p, 2}), ...
                    45, colors(a,:), 'o', 'filled', ...
                    'MarkerFaceAlpha', 1.0, ...
                    'MarkerEdgeColor', [0.2 0.2 0.2], ...
                    'LineWidth', 0.5, ...
                    'DisplayName', labels{a});
            else
                scatter(all_data(:, proj_configs{p, 1}), all_data(:, proj_configs{p, 2}), ...
                    20, colors(a,:), 'o', 'filled', ...
                    'MarkerFaceAlpha', 0.08, ...
                    'MarkerEdgeColor', 'none', ...
                    'HandleVisibility', 'on', ...
                    'DisplayName', labels{a});
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

saveas(fig4, fullfile(output_dir, 'comparison_projections_2d.png'));
saveas(fig4, fullfile(output_dir, 'comparison_projections_2d.fig'));
close(fig4);

fprintf('正在绘制收敛速度与成功率图表...\n');

num_algs = length(algorithms);
n_runs_actual = length(results.(algorithms{1}).runtimes);
all_runtimes = zeros(n_runs_actual, num_algs);
all_success = zeros(n_runs_actual, num_algs);

for a = 1:num_algs
    alg_name = algorithms{a};
    if isfield(results.(alg_name), 'runtimes')
        all_runtimes(:, a) = results.(alg_name).runtimes(:);
        all_success(:, a) = results.(alg_name).success_rates(:) * 100;
    else
        warning('未找到 %s 的 runtime 或 success_rate 数据，请检查 mat 文件。', alg_name);
    end
end

fig_runtime = figure('Position', [100, 100, 600, 450]);
set(gcf, 'Color', 'w');
hold on; grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 12, 'LineWidth', 1.2);

total_iterations = 300;
mean_rt = mean(all_runtimes, 1) / total_iterations;
std_rt = std(all_runtimes, 0, 1) / total_iterations;

for a = 1:num_algs
    bar(a, mean_rt(a), 0.6, 'FaceColor', colors(a,:), 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.8);
end
errorbar(1:num_algs, mean_rt, std_rt, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);

set(gca, 'XTick', 1:num_algs, 'XTickLabel', labels);
ylabel('Average Runtime per Iteration (Seconds)', 'FontWeight', 'bold', 'FontSize', 12);
title('Algorithm Runtime Comparison', 'FontWeight', 'bold', 'FontSize', 14);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
box on;

exportgraphics(fig_runtime, fullfile(output_dir, 'comparison_runtime_bar.png'), 'Resolution', 300);
saveas(fig_runtime, fullfile(output_dir, 'comparison_runtime_bar.fig'));
close(fig_runtime);

fig_success = figure('Position', [150, 150, 600, 450]);
set(gcf, 'Color', 'w');

b2 = boxplot(all_success, 'Labels', labels, 'Colors', 'k');
set(b2, 'LineWidth', 1.5);
h = findobj(gca, 'Tag', 'Box');
for j = 1:length(h)
    patch(get(h(j), 'XData'), get(h(j), 'YData'), colors(num_algs-j+1,:), 'FaceAlpha', 0.6);
end

ylabel('Task Execution Success Rate (%)', 'FontWeight', 'bold', 'FontSize', 12);
title('Task Success Rate Distribution (30 Independent Runs)', 'FontWeight', 'bold', 'FontSize', 14);
grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 12, 'LineWidth', 1.2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
box on;

exportgraphics(fig_success, fullfile(output_dir, 'comparison_success_rate_boxplot.png'), 'Resolution', 300);
saveas(fig_success, fullfile(output_dir, 'comparison_success_rate_boxplot.fig'));
close(fig_success);

fprintf('Runtime and success rate charts saved to %s\n', output_dir);

%% ==== 新增：IGD 与 Spread 的核心指标可视化 ====
fprintf('正在绘制 IGD 与 Spread 箱线图...\n');

all_igd = zeros(n_runs_actual, num_algs);
all_spread = zeros(n_runs_actual, num_algs);

for a = 1:num_algs
    alg_name = algorithms{a};
    if isfield(results.(alg_name), 'igd_values') && isfield(results.(alg_name), 'spread_values')
        all_igd(:, a) = results.(alg_name).igd_values(:);
        all_spread(:, a) = results.(alg_name).spread_values(:);
    else
        warning('未找到 %s 的 IGD 或 Spread 数据！', alg_name);
    end
end

fig_igd = figure('Position', [200, 200, 600, 450]);
set(gcf, 'Color', 'w');
b_igd = boxplot(all_igd, 'Labels', labels, 'Colors', 'k');
set(b_igd, 'LineWidth', 1.5);
h_igd = findobj(gca,'Tag','Box');
for j = 1:length(h_igd)
    patch(get(h_igd(j),'XData'), get(h_igd(j),'YData'), colors(num_algs-j+1,:), 'FaceAlpha', 0.6);
end
ylabel('Inverted Generational Distance (IGD) \downarrow', 'FontWeight', 'bold', 'FontSize', 12);
title('Convergence Evaluation: IGD Metric (Lower is Better)', 'FontWeight', 'bold', 'FontSize', 14);
grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 11, 'LineWidth', 1.2);
exportgraphics(fig_igd, fullfile(output_dir, 'comparison_igd_boxplot.png'), 'Resolution', 300);
saveas(fig_igd, fullfile(output_dir, 'comparison_igd_boxplot.fig'));
close(fig_igd);

fig_spread = figure('Position', [250, 250, 600, 450]);
set(gcf, 'Color', 'w');
b_spread = boxplot(all_spread, 'Labels', labels, 'Colors', 'k');
set(b_spread, 'LineWidth', 1.5);
h_spread = findobj(gca,'Tag','Box');
for j = 1:length(h_spread)
    patch(get(h_spread(j),'XData'), get(h_spread(j),'YData'), colors(num_algs-j+1,:), 'FaceAlpha', 0.6);
end
ylabel('Spread Metric \uparrow', 'FontWeight', 'bold', 'FontSize', 12);
title('Diversity Evaluation: Spread Metric (Higher is Better)', 'FontWeight', 'bold', 'FontSize', 14);
grid on; set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6, 'FontSize', 11, 'LineWidth', 1.2);
exportgraphics(fig_spread, fullfile(output_dir, 'comparison_spread_boxplot.png'), 'Resolution', 300);
saveas(fig_spread, fullfile(output_dir, 'comparison_spread_boxplot.fig'));
close(fig_spread);

fprintf('IGD and Spread charts saved to %s\n', output_dir);
fprintf('  - comparison_igd_boxplot.fig/png (收敛逼近度)\n');
fprintf('  - comparison_spread_boxplot.fig/png (解集分布广度)\n');
end