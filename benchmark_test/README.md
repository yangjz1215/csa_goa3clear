# CEC 2017 Benchmark 测试指南

## 第一步：下载 CEC2017 官方数据包

1. 访问 CEC 2017 官网：
   http://web.myspn.com.tw/cec2017/

2. 下载以下文件：
   - `cec2017_data.tar.gz` (包含所有测试函数的数据文件)
   - `cec17_func.cpp` (源代码)

3. 解压并整理目录结构：
```
benchmark_test/
├── cec2017_data/
│   ├── cec2017/
│   │   ├── input_data/      # Shift/Rotation 数据
│   │   ├── cec17_func.cpp   # 源代码
│   │   └── ...
│   └── input_data/           # 也可能直接放在这里
├── cec17_func.mexw64        # Windows 编译产物
└── ...
```

## 第二步：编译 MEX 接口

### Windows (MATLAB)
```matlab
cd d:\csa_goa3\benchmark_test\cec2017_data\cec2017
mex cec17_func.cpp
```

### Linux 服务器
```bash
cd ~/csa_goa3/benchmark_test/cec2017_data/cec2017
mex cec17_func.cpp
```

编译成功后会在当前目录生成：
- Windows: `cec17_func.mexw64`
- Linux: `cec17_func.mexa64`

## 第三步：快速测试

在 MATLAB 中运行：

```matlab
cd d:\csa_goa3\benchmark_test

% 测试单个函数
[best_fit, best_x, curve] = main_cec(1, 30, 300000, -100, 100);

% 检查收敛曲线
plot(curve);
title('F1 Convergence - cSA-GOA');
xlabel('FES');
ylabel('Fitness (Minimizing)');
```

## 第四步：运行完整实验

```matlab
cd d:\csa_goa3\benchmark_test

% 运行 F1-F30 所有函数
results = run_cec_experiments(1:30, 30, 300000, 30);
```

## 注意事项

1. **数据文件路径**：CEC2017 需要读取 `input_data/` 目录下的数据文件
   - 确保 MATLAB 当前目录在 `cec2017_data/` 或其父目录

2. **维度匹配**：确保 `Dim` 是偶数（因为 N_UAV = Dim/2）

3. **最优值参考**（CEC2017）：
   - F1: 100
   - F2: 1100
   - F3: 700
   - ... (详见 CEC2017 官方文档)

4. **误差计算**：Final Error = f(x) - f(x*)

## 快速故障排除

### Error: mex failed
确保安装了 C++ 编译器：
```matlab
mex -setup C++
```

### Error: input_data not found
确保在正确的目录下运行，或修改 `cec17_func.cpp` 中的路径

### Error: dimension mismatch
确保 `Dim` 是偶数，且与 `cec17_func.cpp` 编译时设置一致