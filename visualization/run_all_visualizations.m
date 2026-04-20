function run_all_visualizations(output_dir)
    fprintf('========== Running All Visualizations ==========\n\n');

    fprintf('1. Generating Ablation Charts...\n');
    try
        plot_ablation_all();
        fprintf('   -> Done\n\n');
    catch ME
        fprintf('   -> Error: %s\n\n', ME.message);
    end

    fprintf('2. Generating Sensitivity Analysis Charts...\n');
    try
        plot_sensitivity();
        fprintf('   -> Done\n\n');
    catch ME
        fprintf('   -> Error (or no data): %s\n\n', ME.message);
    end

    fprintf('3. Generating Comparison Charts...\n');
    try
        plot_comparison_all();
        fprintf('   -> Done\n\n');
    catch ME
        fprintf('   -> Error (or no data): %s\n\n', ME.message);
    end

    fprintf('4. Generating Deployment Topology...\n');
    try
        plot_deployment_topology();
        fprintf('   -> Done\n\n');
    catch ME
        fprintf('   -> Error: %s\n\n', ME.message);
    end

    fprintf('========== All Visualizations Complete ==========\n');
    fprintf('Please check the following directories for your figures:\n');
    fprintf('  - figures/ablation\n');
    fprintf('  - figures/sensitivity\n');
    fprintf('  - figures/comparison\n');
    fprintf('  - figures (Topology)\n');
end