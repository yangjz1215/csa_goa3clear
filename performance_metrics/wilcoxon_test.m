function [p_value, h_stat] = wilcoxon_test(data1, data2)
    [p_value, h_stat] = ranksum(data1, data2);
end

function [p_value, h_stat] = wilcoxon_signed_rank_test(data1, data2)
    diff = data1 - data2;
    diff = diff(diff ~= 0);
    n = length(diff);

    if n == 0
        p_value = 1;
        h_stat = 0;
        return;
    end

    [abs_diff, idx] = sort(abs(diff));
    ranks = tie_correction(abs_diff);

    signed_ranks = zeros(n, 1);
    for i = 1:n
        if diff(idx(i)) > 0
            signed_ranks(i) = ranks(i);
        else
            signed_ranks(i) = -ranks(i);
        end
    end

    W = sum(signed_ranks);

    [~, p_value] = ztest(W, 0, sqrt(n*(n+1)*(2*n+1)/6));

    h_stat = (p_value < 0.05);
end

function ranks = tie_correction(abs_diff)
    n = length(abs_diff);
    ranks = zeros(n, 1);
    i = 1;

    while i <= n
        j = i;
        while j < n && abs_diff(j+1) == abs_diff(i)
            j = j + 1;
        end

        rank_sum = sum(i:j);
        avg_rank = rank_sum / (j - i + 1);
        ranks(i:j) = avg_rank;

        i = j + 1;
    end
end

function [p_value, h_stat] = statistical_summary(data, alpha)
    if nargin < 2
        alpha = 0.05;
    end

    n = length(data);
    mean_val = mean(data);
    std_val = std(data);
    median_val = median(data);

    fprintf('样本数: %d\n', n);
    fprintf('均值: %.4f\n', mean_val);
    fprintf('标准差: %.4f\n', std_val);
    fprintf('中位数: %.4f\n', median_val);
    fprintf('最小值: %.4f\n', min(data));
    fprintf('最大值: %.4f\n', max(data));

    if nargout == 0
        fprintf('95%%置信区间: [%.4f, %.4f]\n', ...
            mean_val - 1.96*std_val/sqrt(n), ...
            mean_val + 1.96*std_val/sqrt(n));
    end
end
