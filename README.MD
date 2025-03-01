# Deformation Tranfer 算法

Deformation Transfer 算法的朴素实现，基于 Eigen 和 Openmesh。实现只包括形变迁移损失部分，这意味着输入的源模型和目标模型需要是同拓扑的（相同的点和面）。

使用方法：

1. 配置 Eigen 库和 Openmesh 库。

2. 构建项目，运行

```shell
./dt filename_s0 filename_s1 filename_t0 filename_t1
```

其中， 前两个文件对应源模型的动作A, B，第三个文件对应目标模型的动作A，第四个文件对应结果（目标模型的动作B）输出的路径。

或将模型文件保存在 `./objects/s0( or s1, t0).obj`，直接运行 `./run.sh` 即可。
