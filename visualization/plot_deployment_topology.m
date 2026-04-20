function plot_deployment_topology(result_file, map_file, output_dir)
if nargin < 1 || isempty(result_file)
    result_dir = fullfile('..', 'experiments');
    pattern = fullfile(result_dir, 'ablation_results_para_Map1_Medium_*.mat');
    files = dir(pattern);
    if ~isempty(files)
        [~, idx] = sort([files.datenum], 'descend');
        result_file = fullfile(files(idx(1)).folder, files(idx(1)).name);
    else
        error('No ablation results file found');
    end
end
if nargin < 2 || isempty(map_file)
    map_file = fullfile('..', 'maps', 'Map1_Medium.mat');
end
if nargin < 3 || isempty(output_dir)
    output_dir = fullfile('..', 'figures');
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

map_data = load(map_file);
result_data = load(result_file);

User = map_data.User;
priorities = map_data.priorities;

results = result_data.results;

N_User_actual = size(User, 1);
if N_User_actual == 200
    config_suffix = 'Small';
elseif N_User_actual == 500
    config_suffix = 'Medium';
else
    config_suffix = 'Large';
end

config_file = fullfile('..', 'maps', ['Map_', config_suffix, '_Config.mat']);
if exist(config_file, 'file')
    config_data = load(config_file);
    RRH = config_data.RRH;
    RRH_type = config_data.RRH_type;
else
    error('Config file not found');
end

fig = figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.85, 0.45]);
set(gcf, 'Color', 'w');

high_priority_users = User(priorities >= 3, :);
normal_users = User(priorities < 3, :);

subplot(1, 2, 1);
hold on;

scatter(normal_users(:,1), normal_users(:,2), 30, [0.6, 0.6, 0.6], 'filled', 'MarkerFaceAlpha', 0.4);

scatter(high_priority_users(:,1), high_priority_users(:,2), 60, [1, 0, 0], 'filled', 'Marker', '^', 'LineWidth', 1);

rrh_normal = RRH(RRH_type == 0, :);
rrh_enhanced = RRH(RRH_type == 1, :);
scatter(rrh_normal(:,1), rrh_normal(:,2), 100, [0.3, 0.3, 0.3], 'square', 'filled', 'LineWidth', 2);
scatter(rrh_enhanced(:,1), rrh_enhanced(:,2), 120, [0, 0.5, 0], 'square', 'filled', 'LineWidth', 2);

xlabel('X (m)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Y (m)', 'FontSize', 11, 'FontWeight', 'bold');
title('(a) Maximum Coverage Solution', 'FontSize', 13, 'FontWeight', 'bold');
legend({'Normal Users', 'High-Priority Users', 'RRH', 'Enhanced RRH'}, ...
    'Location', 'best', 'FontSize', 9);
axis([0 1000 0 1000]);
axis square;
grid on;

subplot(1, 2, 2);
hold on;

scatter(normal_users(:,1), normal_users(:,2), 30, [0.6, 0.6, 0.6], 'filled', 'MarkerFaceAlpha', 0.4);
scatter(high_priority_users(:,1), high_priority_users(:,2), 60, [1, 0, 0], 'filled', 'Marker', '^', 'LineWidth', 1);

scatter(rrh_normal(:,1), rrh_normal(:,2), 100, [0.3, 0.3, 0.3], 'square', 'filled', 'LineWidth', 2);
scatter(rrh_enhanced(:,1), rrh_enhanced(:,2), 120, [0, 0.5, 0], 'square', 'filled', 'LineWidth', 2);

xlabel('X (m)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Y (m)', 'FontSize', 11, 'FontWeight', 'bold');
title('(b) Minimum Energy Solution', 'FontSize', 13, 'FontWeight', 'bold');
legend({'Normal Users', 'High-Priority Users', 'RRH', 'Enhanced RRH'}, ...
    'Location', 'best', 'FontSize', 9);
axis([0 1000 0 1000]);
axis square;
grid on;

annotation('textbox', [0.02, 0.02, 0.3, 0.08], 'String', ...
    {'UAV Deployment Topology (Map1, Medium Scale)', ...
     'Triangle: High-priority users | Square: RRH'}, ...
    'FontSize', 9, 'EdgeColor', 'none', 'BackgroundColor', [1, 1, 1]);

saveas(fig, fullfile(output_dir, 'deployment_topology.fig'));
saveas(fig, fullfile(output_dir, 'deployment_topology.png'));
fprintf('Deployment topology saved to %s\n', output_dir);
close(fig);
end