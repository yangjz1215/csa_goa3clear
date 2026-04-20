function setupAndGenerateMaps()
    fprintf('========== 初始化并生成固定测试地图 ==========\n\n');

    project_dir = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(project_dir));

    cd(fullfile(project_dir, 'maps'));
    fprintf('当前目录: %s\n', pwd);

    generateMaps();

    fprintf('\n========== 地图生成完成 ==========\n');
    fprintf('现在可以在实验中使用这些固定地图了。\n');
end