function generateMaps()
    fprintf('========== 生成固定测试地图 ==========\n');

    Lb = [0, 0];
    Ub = [1000, 1000];
    center = [500, 500];

    scenario_configs = {
        'Small', 200, 8;
        'Medium', 500, 15;
        'Large', 1000, 25
    };

    map_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'maps');
    if ~exist(map_dir, 'dir')
        mkdir(map_dir);
    end

    for s_idx = 1:size(scenario_configs, 1)
        scenario_name = scenario_configs{s_idx, 1};
        N_User = scenario_configs{s_idx, 2};
        N_UAV = scenario_configs{s_idx, 3};

        fprintf('\n--- 生成场景: %s (%d用户, %dUAV) ---\n', scenario_name, N_User, N_UAV);

        ratios = [0.25, 0.25, 0.25, 0.25];
        counts = floor(ratios * N_User);
        remainder = N_User - sum(counts);
        extra_indices = randperm(length(ratios), remainder);
        for i = 1:remainder
            counts(extra_indices(i)) = counts(extra_indices(i)) + 1;
        end
        task_levels = [];
        for i = 1:length(ratios)
            task_levels = [task_levels; i * ones(counts(i), 1)];
        end
        priorities = task_levels(randperm(N_User));

        N_RRH = 10;
        N_eRRH = 4;
        RRH = Lb + (Ub - Lb) .* rand(N_RRH, 2);
        RRH_type = zeros(N_RRH, 1);
        RRH_type(1:N_eRRH) = 1;

        N_eUAV = floor(N_UAV * 0.3);
        UAV_type = zeros(N_UAV, 1);
        UAV_type(1:N_eUAV) = 1;

        D_max = 500; D_min = 100;
        D = ((D_max-D_min)*rand(N_User,1)+D_min)*8192;
        C_max = 1.0; C_min = 0.5;
        C = ((C_max-C_min)*rand(N_User,1)+C_min)*10^9;
        DT_max = 1.2; DT_min = 0.8;
        DT = ((DT_max-DT_min)*rand(N_User,1)+DT_min);

        save(fullfile(map_dir, sprintf('Map_%s_Config.mat', scenario_name)), ...
            'N_User', 'N_UAV', 'N_RRH', 'N_eRRH', 'N_eUAV', ...
            'RRH', 'RRH_type', 'UAV_type', 'priorities', ...
            'D', 'C', 'DT', 'Lb', 'Ub', 'center');

        fprintf('  保存配置: Map_%s_Config.mat\n', scenario_name);
    end

    map_types = {
        1, 'Map1', 'Uniform', '均匀分布';
        2, 'Map2', 'BiCenter', '双中心聚集';
        3, 'Map3', 'EdgeScatter', '边缘散落';
        4, 'Map4', 'SingleCenter', '单中心密集';
        5, 'Map5', 'MultiCenter', '多中心分散'
    };

    for m_idx = 1:size(map_types, 1)
        map_num = map_types{m_idx, 1};
        map_prefix = map_types{m_idx, 2};
        generator_name = map_types{m_idx, 3};
        map_desc = map_types{m_idx, 4};

        fprintf('\n--- 生成地图类型: %s (%s) ---\n', map_prefix, map_desc);

        for s_idx = 1:size(scenario_configs, 1)
            scenario_name = scenario_configs{s_idx, 1};
            N_User = scenario_configs{s_idx, 2};

            switch generator_name
                case 'Uniform'
                    User = Lb + (Ub - Lb) .* rand(N_User, 2);
                case 'BiCenter'
                    User = generateBiCenter(N_User, center, [200, 100]);
                case 'EdgeScatter'
                    User = generateEdgeScatter(N_User, Lb, Ub, center);
                case 'SingleCenter'
                    User = center + 150 * randn(N_User, 2);
                    User = max(0, min(1000, User));
                case 'MultiCenter'
                    User = generateMultiCenter(N_User, center, 4);
            end

            counts = floor([0.25, 0.25, 0.25, 0.25] * N_User);
            remainder = N_User - sum(counts);
            extra_indices = randperm(4, remainder);
            for i = 1:remainder
                counts(extra_indices(i)) = counts(extra_indices(i)) + 1;
            end
            task_levels = [];
            for i = 1:4
                task_levels = [task_levels; i * ones(counts(i), 1)];
            end
            priorities = task_levels(randperm(N_User));

            filename = sprintf('%s_%s.mat', map_prefix, scenario_name);
            save(fullfile(map_dir, filename), 'User', 'priorities', 'N_User', 'Lb', 'Ub');

            fprintf('  [%s] %s: %s (N_User=%d)\n', scenario_name, map_prefix, filename, N_User);
        end
    end

    fprintf('\n========== 地图生成完成 ==========\n');
    fprintf('地图保存位置: %s\n', map_dir);
    fprintf('\n生成的地图文件:\n');
    dir_files = dir(fullfile(map_dir, 'Map*.mat'));
    for i = 1:length(dir_files)
        fprintf('  - %s\n', dir_files(i).name);
    end
end

function User = generateBiCenter(N_User, center, radii)
    User = zeros(N_User, 2);
    half = floor(N_User / 2);
    for i = 1:half
        r = radii(1) * sqrt(rand());
        theta = 2 * pi * rand();
        User(i, :) = center + r * [cos(theta), sin(theta)];
    end
    offset = [200, 0];
    for i = half+1:N_User
        r = radii(2) * sqrt(rand());
        theta = 2 * pi * rand();
        User(i, :) = center + offset + r * [cos(theta), sin(theta)];
    end
    if mod(N_User, 2) == 1
        User(end, :) = center + (rand(1,2) - 0.5) * 100;
    end
    User = max(0, min(1000, User));
end

function User = generateEdgeScatter(N_User, Lb, Ub, center)
    User = zeros(N_User, 2);
    margin = 150;
    n_edge = floor(N_User * 0.7);
    for i = 1:n_edge
        edge = randi(4);
        switch edge
            case 1
                User(i, :) = [Lb(1) + margin + rand() * (Ub(1) - 2*margin), Lb(2) + rand() * (Ub(2) - 2*margin)];
            case 2
                User(i, :) = [Lb(1) + rand() * (Ub(1) - 2*margin), Lb(2) + margin + rand() * (Ub(2) - 2*margin)];
            case 3
                User(i, :) = [Ub(1) - margin - rand() * (Ub(1) - 2*margin), Lb(2) + rand() * (Ub(2) - 2*margin)];
            case 4
                User(i, :) = [Lb(1) + rand() * (Ub(1) - 2*margin), Ub(2) - margin - rand() * (Ub(2) - 2*margin)];
        end
    end
    n_corner = N_User - n_edge;
    corner_size = 100;
    corners = [
        Lb + [0, 0];
        Lb + [Ub(1)-corner_size, 0];
        Lb + [0, Ub(2)-corner_size];
        Ub - [corner_size, corner_size]
    ];
    for i = 1:n_corner
        corner_idx = randi(4);
        User(n_edge + i, :) = corners(corner_idx, :) + rand(1,2) * corner_size;
    end
    User = max(0, min(1000, User));
end

function User = generateMultiCenter(N_User, center, n_centers)
    User = zeros(N_User, 2);
    spacing = 300;
    offsets = zeros(n_centers, 2);
    for c = 1:n_centers
        row = floor((c-1) / 2);
        col = mod(c-1, 2);
        offsets(c, :) = (col - 0.5) * spacing + (row - 0.5) * spacing;
    end
    per_center = floor(N_User / n_centers);
    idx = 1;
    for c = 1:n_centers
        start_idx = (c-1) * per_center + 1;
        end_idx = c * per_center;
        if c == n_centers
            end_idx = N_User;
        end
        cluster_size = length(start_idx:end_idx);
        cluster_center = center + offsets(c, :);
        User(start_idx:end_idx, :) = cluster_center + 80 * randn(cluster_size, 2);
    end
    User = max(0, min(1000, User));
end