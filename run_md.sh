#!/usr/bin/env bash
# 自动化流程：从 PDB 拆分、拓扑生成，到 GROMACS MD 模拟和 MMPBSA 能量计算
set -e  # 一出错就停止执行

# 判断是否传入参数
if [ -z "$1" ]; then
  echo "❌ 错误：未提供输入 PDB 文件。请这样调用：bash run.sh ER_ZEN.pdb"
  exit 1
fi
# 判断是否传入参数
if [ -z "$1" ]; then
  echo "❌ 错误：未提供输入 PDB 文件。"
  echo "用法: bash run.sh <input.pdb> <gmx_command>"
  exit 1
fi

if [ -z "$2" ]; then
  echo "❌ 错误：未提供 GMX 执行路径。"
  echo "用法: bash run.sh <input.pdb> <gmx_command>"
  exit 1
fi
INPUT_PDB="$1"
GMX="$2"

echo "==> 1. 拆分 PDB 为蛋白与配体"
python split_pdb.py "$INPUT_PDB"  # 输出: protein.pdb, ligand.pdb

echo "==> 2. 加载 Intel/GROMACS 环境"
#source /opt/intel/env.sh
#GMX="/opt/software/gromacs2021-GPU/bin/gmx_mpi"

echo "==> 3. 生成蛋白拓扑"
$GMX pdb2gmx -f protein.pdb -o processed.gro -water spc <<EOF
6
EOF

echo "==> 4. 生成配体拓扑（AmberTools + ACPYPE）"

antechamber -i ligand.pdb -fi pdb -o ligand.mol2 -fo mol2 -c bcc -s 2 -nc 0
parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod
acpype -i ligand.mol2 -b Ligand

# 清理中间文件
mv Ligand.acpype/Ligand_GMX.gro .
mv Ligand.acpype/Ligand_GMX.itp .
rm -rf Ligand.acpype ligand.frcmod ligand.mol2 sqm* ANTECHAMBER* ATOMTYPE.INF

# 删除 gro 中非标准行，仅保留原子坐标
sed '1,2d' Ligand_GMX.gro | sed '$d' > Ligand_GMX_1.gro

echo "==> 5. 合并蛋白与配体坐标"
COMPLEX=processed.gro
LIGAND=Ligand_GMX_1.gro
OUT=complex.gro
box_line=$(tail -n1 "$COMPLEX")

head -n -1 "$COMPLEX" > tmp1.gro
cat "$LIGAND" >> tmp1.gro
echo "$box_line" >> tmp1.gro
sed -i '/^\s*$/d' tmp1.gro

title=$(head -n1 tmp1.gro)
atom_lines=$(( $(wc -l < tmp1.gro) - 3 ))
{
  echo "$title"
  printf "%5d\n" "$atom_lines"
  tail -n+3 tmp1.gro
} > "$OUT"

rm tmp1.gro Ligand_GMX_1.gro ligand.pdb processed.gro protein.pdb 
echo "  ✅ 已生成 $OUT （原子数: $atom_lines）"

echo "==> 6. 更新 topol.top 以包含配体"
TOPOL=topol.top
cp "$TOPOL" "$TOPOL.bak"

awk '
  /#include "amber99sb-ildn.ff\/forcefield.itp"/ && !done {
    print
    print "; Include ligand topology"
    print "#include \"Ligand_GMX.itp\""
    done=1
    next
  }
  { print }
' "$TOPOL.bak" > "$TOPOL"

# 追加分子计数
grep -q "^Ligand" "$TOPOL" || echo "Ligand                  1" >> "$TOPOL"
echo "  ✅ 已更新 $TOPOL"

echo "==> 7. 创建盒子、居中复合物、加水"
$GMX editconf -f "$OUT" -o complex_box.gro -bt dodecahedron -d 2.0
python center_gro_precisely.py complex_box.gro
$GMX solvate -cp complex_box_centered.gro -cs spc216.gro -p topol.top -o solv.gro
rm complex_box_centered.gro complex_box.gro 

echo "==> 8. 加离子并能量最小化"
$GMX grompp -f em.mdp -c solv.gro -p topol.top -o em.tpr -maxwarn 1
$GMX genion -s em.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral  <<EOF
15
EOF
$GMX grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
$GMX mdrun -deffnm em

echo "==> 9. 添加配体位置约束"
$GMX genrestr -f Ligand_GMX.gro -o posre_LI.itp -fc 1000 1000 1000 <<EOF
0
EOF
sed -i '/LI_GMX\.itp/a\
; Ligand position restraints\n#ifdef POSRES\n#include \"posre_Ligand.itp\"\n#endif' topol.top

echo "==> 10. 构建 index 文件（合并蛋白+配体）"
$GMX make_ndx -f em.gro -o index.ndx <<EOF
1 | 19
q
EOF

echo "==> 11. NVT 平衡"
$GMX grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
$GMX mdrun -deffnm nvt -v 

echo "==> 12. NPT 平衡"
$GMX grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr
$GMX mdrun -deffnm npt -v 


echo "==> 13. 生产模拟 50 ns"
$GMX grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_1.tpr
$GMX mdrun -deffnm md_0_1 -v -nb gpu -pme gpu -bonded gpu
rm npt.cpt npt.edr npt.gro npt.log npt.tpr npt.trr nvt.edr
rm em.tpr.1 complex.gro em.edr em.gro em.log em.tpr em.trr topol.top.bak
rm nvt.trr nvt.tpr nvt.log nvt.cpt nvt_prev.cpt nvt.gro solv.gro solv_ions.gro 
#rm "#em.tpr.1#" "#topol.top.1#" "#topol.top.2#"
