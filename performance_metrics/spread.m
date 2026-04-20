function spread_value = spread(obtained_front, true_front)
    if size(obtained_front, 1) < 2
        spread_value = 0;
        return;
    end

    n_points = size(obtained_front, 1);
    n_objectives = size(obtained_front, 2);

    if nargin < 2 || isempty(true_front)
        [sorted_obtained, ~] = sortrows(obtained_front, 1);
        extreme_first = sorted_obtained(1, :);
        extreme_last = sorted_obtained(end, :);
    else
        [sorted_true, ~] = sortrows(true_front, 1);
        extreme_first = sorted_true(1, :);
        extreme_last = sorted_true(end, :);
    end

    [sorted_obtained, idx_sort] = sortrows(obtained_front, 1);

    d_extreme_first = inf;
    for i = 1:n_points
        dist = norm(sorted_obtained(i, :) - extreme_first);
        if dist < d_extreme_first
            d_extreme_first = dist;
        end
    end

    d_extreme_last = inf;
    for i = 1:n_points
        dist = norm(sorted_obtained(i, :) - extreme_last);
        if dist < d_extreme_last
            d_extreme_last = dist;
        end
    end

    d_i = zeros(n_points - 1, 1);
    for i = 1:n_points - 1
        d_i(i) = norm(sorted_obtained(i, :) - sorted_obtained(i + 1, :));
    end

    d_avg = mean(d_i);

    numerator = d_extreme_first + d_extreme_last + sum(abs(d_i - d_avg));
    denominator = d_extreme_first + d_extreme_last + (n_points - 1) * d_avg;

    if denominator == 0
        spread_value = 0;
    else
        spread_value = numerator / denominator;
    end
end
