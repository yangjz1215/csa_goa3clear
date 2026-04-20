function metrics = calculate_all_metrics(obtained_front, true_front, reference_point)
    if nargin < 3
        reference_point = [];
    end
    if nargin < 2
        true_front = [];
    end

    metrics = struct();

    if ~isempty(obtained_front)
        if isempty(reference_point)
            if isempty(true_front)
                reference_point = max(obtained_front, [], 1) * 1.1;
            else
                reference_point = max([obtained_front; true_front], [], 1) * 1.1;
            end
        end
        metrics.hv = hypervolume(obtained_front, reference_point);
    else
        metrics.hv = 0;
    end

    if ~isempty(obtained_front) && ~isempty(true_front)
        max_objs = max(true_front, [], 1);
        min_objs = min(true_front, [], 1);
        range_objs = max_objs - min_objs + 1e-6;

        norm_true_front = true_front;
        norm_algo_front = obtained_front;

        norm_true_front(:, 1) = (max_objs(1) - true_front(:, 1)) / range_objs(1);
        norm_algo_front(:, 1) = (max_objs(1) - obtained_front(:, 1)) / range_objs(1);

        norm_true_front(:, 2:3) = (true_front(:, 2:3) - min_objs(2:3)) ./ range_objs(2:3);
        norm_algo_front(:, 2:3) = (obtained_front(:, 2:3) - min_objs(2:3)) ./ range_objs(2:3);

        metrics.igd = igd(norm_algo_front, norm_true_front);
        metrics.gd = gd(norm_algo_front, norm_true_front);
    else
        metrics.igd = nan;
        metrics.gd = nan;
        norm_algo_front = [];
        norm_true_front = [];
    end

    if ~isempty(obtained_front) && ~isempty(norm_algo_front)
        metrics.spread = spread(norm_algo_front, norm_true_front);
    else
        metrics.spread = nan;
    end

    if ~isempty(obtained_front)
        metrics.n_solutions = size(obtained_front, 1);
    else
        metrics.n_solutions = 0;
    end
end
