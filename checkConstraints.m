function X_checked = checkConstraints(X, params)
    N = params.N_UAV;
    dim = N * 3;
    
    Ub_rep = repmat(params.Ub, 1, N);
    Lb_rep = repmat(params.Lb, 1, N);
    X_checked = max(min(X, Ub_rep), Lb_rep);
    
    if isfield(params, 'd_safe')
        pos = reshape(X_checked, 3, N)';
        for i = 1:N
            for j = i+1:N
                if norm(pos(i,:) - pos(j,:)) < params.d_safe
                    pos(i,:) = pos(i,:) + (rand(1,3)-0.5) * params.d_safe;
                end
            end
        end
        X_checked = reshape(pos', 1, dim);
        X_checked = max(min(X_checked, Ub_rep), Lb_rep);
    end
end
