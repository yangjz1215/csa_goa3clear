function fitness = benchmark_wrapper(X, func_num)
    x_col = X(:);
    fitness = cec17_func(x_col, func_num);
end