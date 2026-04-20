function hv = hypervolume(objectives, reference_point)
    if nargin < 2
        reference_point = max(objectives, [], 1) * 1.1;
    end

    n_points = size(objectives, 1);
    n_objectives = size(objectives, 2);

    if n_points == 0
        hv = 0;
        return;
    end

    obj_norm = normalize(objectives, reference_point);

    if n_objectives == 2
        hv = hv_2d(obj_norm);
    elseif n_objectives == 3
        hv = hv_3d(obj_norm);
    else
        hv = hv_nd(obj_norm, n_objectives);
    end
end

function obj_norm = normalize(objectives, reference_point)
    obj_norm = objectives;
    for m = 1:size(objectives, 2)
        obj_norm(:, m) = reference_point(m) - objectives(:, m);
    end
end

function hv = hv_2d(objectives)
    [~, idx] = sort(objectives(:, 1), 'descend');
    objectives = objectives(idx, :);

    hv = 0;
    prev_y = 0;

    for i = 1:size(objectives, 1)
        width = objectives(i, 1);
        height = objectives(i, 2) - prev_y;
        hv = hv + width * height;
        prev_y = objectives(i, 2);
    end
end

function hv = hv_3d(objectives)
    [~, idx] = sortrows(objectives, [-1 -2 -3]);
    objectives = objectives(idx, :);

    hv = 0;
    included = false(size(objectives, 1), 1);

    for i = 1:size(objectives, 1)
        included(i) = true;
        current_points = objectives(included, :);

        [~, idx2] = sortrows(current_points, [-2 -3]);
        current_points = current_points(idx2, :);

        contribution = current_points(1, 1);

        if size(current_points, 1) > 1
            area_2d = 0;
            prev_z = 0;
            for j = 1:size(current_points, 1)
                width = current_points(j, 2);
                height = current_points(j, 3) - prev_z;
                area_2d = area_2d + width * height;
                prev_z = current_points(j, 3);
            end
            contribution = contribution * area_2d;
        else
            contribution = contribution * current_points(1, 2) * current_points(1, 3);
        end

        hv = hv + contribution;
    end
end

function hv = hv_nd(objectives, n_objectives)
    [~, idx] = sortrows(objectives, -1);
    objectives = objectives(idx, :);

    hv = 0;

    for i = 1:size(objectives, 1)
        if i == 1
            volume = prod(objectives(i, :));
            hv = hv + volume;
        else
            for j = 1:i-1
                diff = objectives(i, :) - objectives(j, :);
                if all(diff <= 0)
                    volume = prod(objectives(i, :));
                    hv = hv + volume;
                    break;
                end
            end
        end
    end
end
