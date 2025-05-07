# 自动化 MD + MMPBSA 流程

本项目提供一键化脚本 `run.sh`，完成从 PDB 拆分、拓扑生成，到 GROMACS 分子动力学（MD）模拟与 MMPBSA 能量计算的全流程。

## 目录结构

```
.
├── run.sh                    # 主脚本
├── split_pdb.py             # 拆分 PDB（蛋白/配体）
├── center_gro_precisely.py  # 精确中心化坐标
├── em.mdp                    # 能量最小化参数
├── nvt.mdp                   # NVT 平衡参数
├── npt.mdp                   # NPT 平衡参数
├── md.mdp                    # 生产模拟参数
├── topol.top                 # GROMACS 拓扑文件（需预先准备好蛋白力场）
└── spc216.gro                # 水分子坐标文件（GROMACS 自带）
```

## 功能简介

1. **拆分 PDB**：将输入的复合物 PDB 拆分为 `protein.pdb`、`ligand.pdb`  
2. **生成蛋白拓扑**：调用 GROMACS `pdb2gmx`  
3. **生成配体拓扑**：通过 AmberTools (`antechamber`、`parmchk2`) + `acpype`  
4. **合并坐标**：合并蛋白/配体坐标，生成 `complex.gro`  
5. **更新拓扑**：在 `topol.top` 中引入配体 `.itp`，并添加分子计数  
6. **装盒、加水、加离子**：`editconf` → `solvate` → `genion` → 能量最小化  
7. **添加配体位置约束**：`genrestr` + 更新 `posre_Ligand.itp`  
8. **生成索引**：`make_ndx` 合并蛋白+配体的 index  
9. **NVT / NPT 平衡**  
10. **生产模拟 50 ns**  
11. **MMPBSA**（可在脚本末尾补充调用）

## 环境与依赖

- **操作系统**：Linux（建议 Ubuntu 18.04+ / CentOS7+）
- **Shell**：bash（`#!/usr/bin/env bash`）
- **Python**：3.6+，需安装 `split_pdb.py`、`center_gro_precisely.py` 所需依赖
- **GROMACS**：2021 或更高，需支持 MPI/GPU（如 `gmx_mpi`）
- **AmberTools**：含 `antechamber`、`parmchk2`
- **ACPYPE**：自动生成 GROMACS 配体拓扑
- **其他**：`sed`、`awk`、`grep` 等常用 GNU 工具

```bash
# 示例：Ubuntu 下快速安装
sudo apt update
sudo apt install -y gromacs python3-pip ambertools
pip3 install acpype      # 若未内置
```

> **注意**：如果使用 Intel OneAPI 或其他定制环境，可在脚本中取消注释：  
> ```bash
> # source /opt/intel/env.sh
> # GMX="/opt/software/gromacs2021-GPU/bin/gmx_mpi"
> ```

## 配置项

| 变量               | 含义                                                        | 示例                       |
|--------------------|-------------------------------------------------------------|----------------------------|
| `INPUT_PDB`        | 输入复合物 PDB 文件                                         | `ER_ZEN.pdb`               |
| `GMX`              | GROMACS 可执行命令或完整路径（含 pdb2gmx、grompp、mdrun）   | `gmx` 或 `/opt/gromacs/bin/gmx_mpi` |
| `em.mdp`, `nvt.mdp` 等 | GROMACS 模拟参数文件，放在当前目录下，文件名不可更改     | —                          |

脚本会依次执行，期间生成的中间文件及目录结构：

```
protein.pdb, ligand.pdb
processed.gro, Ligand_GMX.gro, Ligand_GMX.itp
complex.gro, complex_box.gro, complex_box_centered.gro, solv.gro, solv_ions.gro
em.*, nvt.*, npt.*, md_0_1.*
topol.top（已自动更新 .itp 和分子计数）
index.ndx, posre_Ligand.itp
```

## 使用方法

```bash
bash run.sh <input.pdb> <gmx_command>
```

- **用法示例**：

  ```bash
  bash run.sh ER_ZEN.pdb gmx_mpi
  ```

- **参数说明**：
  1. `<input.pdb>`：复合物结构文件
  2. `<gmx_command>`：GROMACS 命令前缀（如 `gmx`、`gmx_mpi` 或绝对路径）

脚本执行过程中如遇错误会立即停止，并打印 ❌ 错误提示。务必根据提示安装/配置缺失依赖或修正命令。

## 输出结果

- 最终生成 `md_0_1.*`（生产模拟 50 ns）  
- 可在脚本末尾追加 MMPBSA 步骤，例如调用 `g_mmpbsa` 或 `MMPBSA.py`  

## 常见问题

- **未找到 GROMACS 命令**  
  确认 `GMX` 变量正确，或在脚本中指定完整路径。  
- **ACPYPE 生成失败**  
  检查 AmberTools、acpype 是否已安装；配体电荷 (`-nc`) 参数需与分子实际电荷一致。  
- **拓扑冲突**  
  确保 `topol.top` 中 forcefield 路径正确，并且已包含 `Ligand_GMX.itp`。

## 扩展

- **MMPBSA 计算**：可在脚本末尾添加计算命令  
- **参数调优**：`*.mdp` 文件中可根据需求修改时间步长、I/O 频率等
