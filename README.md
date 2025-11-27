# ABS（吸收边界条件）图表复现代码

此仓库包含用于复现相关论文中所有图表的MATLAB代码。

**项目主页**: [CREWES Research Consortium - Free Software](https://www.crewes.org/ResearchLinks/FreeSoftware/)

> 注意：代码中使用了CREWES工具箱中的 `spongeABC` 和 `plotimage` 函数来实现吸收边界条件和绘制图像。请确保您的MATLAB环境已正确配置CREWES工具箱[3](@ref)。

## 图表复现指南

下表汇总了复现每个图表所需运行的具体脚本文件。

| 目标图表 | 所需运行的脚本文件 | 后续整合脚本 |
| :--- | :--- | :--- |
| **Figure 1** | 直接运行 `figure1a.m`, `figure1b.m`, ..., `figure1f.m` | - |
| **Figure 2** | 1. 运行 `figure2_M1_1D.m`, `figure2_M1_1D_38.m`<br>2. 运行 `figure2_M3_BalancedABS_1D.m`, `figure2_M3_BalancedABS_1D_r38.m`, `figure2_M3_Non_BalancedABS_1D.m`, `figure2_M3_Non_BalancedABS_1D_r38.m`<br>3. 运行 `figure2_M4_BalancedABS_1D.m`, `figure2_M4_BalancedABS_1D_r38.m`, `figure2_M4_Non_BalancedABS_1D_r19.m`, `figure2_M4_Non_BalancedABS_1D_r38.m` | 1. 运行 `figure2_ab_compare_plot_2.m` (生成图2a, 2b)<br>2. 运行 `figure2_cd_compare_plot_2.m` (生成图2c, 2d)<br>3. 运行 `figure2_ef_compare_plot.m` (生成图2e, 2f) |
| **Figure 4**<br>(Figure 3 为速度模型) | 运行 `figure3_Tra.m`, `figure3b.m`, `figure3c.m`, `figure3d.m`, `figure3e.m`, `figure3f.m` | 运行 `figure3CompareResult.m` |
| **Figure 6**<br>(Figure 5 为速度模型) | 运行 `figure5a.m`, `figure5b.m`, `figure5c.m`, `figure5d.m`, `figure5e.m`, `figure5f.m` | 运行 `figure5plot.m` |
| **Figure 8**<br>(Figure 7 为速度模型) | 运行 `figure8a.m`, `figure8b.m`, `figure8c.m`, `figure8d.m`, `figure8e.m`, `figure8f.m` | 运行 `figure8CompareResult.m` |
| **Figure 10**<br>(Figure 9 为速度模型) | 运行 `figure10a.m`, `figure10b.m`, `figure10c.m`, `figure10d.m`, `figure10e.m`, `figure10f.m` | 运行 `figure10CompareResult.m` |

## 使用说明

1.  **环境准备**: 确保您已安装MATLAB，并已正确设置CREWES工具箱的路径。
2.  **运行脚本**: 根据上表，按照顺序运行指定脚本。通常需要先运行所有生成数据的脚本，最后运行整合绘图脚本。
3.  **结果查看**: 生成的图表将显示在MATLAB的图形窗口中，并可能自动保存为图像文件（具体取决于脚本设置）。

## 注意事项

-   请严格按照上述顺序运行脚本，因为后续脚本可能依赖于前面脚本生成的数据文件。
-   如果遇到任何错误，请首先检查CREWES工具箱的函数路径是否正确添加至MATLAB的搜索路径中[3](@ref)。
