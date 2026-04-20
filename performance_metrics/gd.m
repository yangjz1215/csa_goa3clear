function gd_value = gd(obtained_front, true_front)
    if isempty(obtained_front)
        gd_value = inf;
        return;
    end

    n_obtained = size(obtained_front, 1);

    min_distances = zeros(n_obtained, 1);

    for i = 1:n_obtained
        obtained_point = obtained_front(i, :);
        min_dist = inf;

        for j = 1:size(true_front, 1)
            true_point = true_front(j, :);
            dist = norm(obtained_point - true_point);
            if dist < min_dist
                min_dist = dist;
            end
        end

        min_distances(i) = min_dist;
    end

    gd_value = sqrt(sum(min_distances.^2)) / n_obtained;
end
