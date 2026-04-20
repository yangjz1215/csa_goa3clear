function [archive, is_updated] = updateParetoArchive3D(archive, new_pos, util, lat, nrg)
    is_updated = false;
    new_obj = [-util, lat, nrg];

    is_dominated = false;
    indices_to_remove = [];

    for i = 1:length(archive)
        old_obj = [-archive(i).Utility, archive(i).Latency, archive(i).Energy];

        if all(old_obj <= new_obj) && any(old_obj < new_obj)
            is_dominated = true;
            break;
        end

        if all(new_obj <= old_obj) && any(new_obj < old_obj)
            indices_to_remove = [indices_to_remove, i];
        end
    end

    if ~is_dominated
        archive(indices_to_remove) = [];
        archive(end+1) = struct('UAV_pos', new_pos, 'Utility', util, 'Latency', lat, 'Energy', nrg);
        is_updated = true;
    end

    MAX_ARCHIVE_SIZE = 300;
    if length(archive) > MAX_ARCHIVE_SIZE
        arch_U = [archive.Utility];
        arch_L = [archive.Latency];
        arch_E = [archive.Energy];

        [~, max_u_idx] = max(arch_U);
        [~, min_l_idx] = min(arch_L);
        [~, min_e_idx] = min(arch_E);
        extreme_indices = unique([max_u_idx, min_l_idx, min_e_idx]);

        all_indices = 1:length(archive);
        other_indices = setdiff(all_indices, extreme_indices);

        num_to_keep = MAX_ARCHIVE_SIZE - length(extreme_indices);
        keep_others = other_indices(randperm(length(other_indices), num_to_keep));

        final_keep_indices = [extreme_indices, keep_others];
        archive = archive(final_keep_indices);
    end
end