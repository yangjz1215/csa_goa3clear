function igd_value = igd(obtained_front, true_front)
    if isempty(obtained_front)
        igd_value = inf;
        return;
    end

    n_true = size(true_front, 1);
    n_obtained = size(obtained_front, 1);

    min_distances = zeros(n_true, 1);

    for i = 1:n_true
        true_point = true_front(i, :);
        min_dist = inf;

        for j = 1:n_obtained
            obtained_point = obtained_front(j, :);
            dist = norm(true_point - obtained_point);
            if dist < min_dist
                min_dist = dist;
            end
        end

        min_distances(i) = min_dist;
    end

    igd_value = mean(min_distances);
end
