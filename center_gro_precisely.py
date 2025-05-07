import numpy as np
import sys
import os

def read_gro(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    title = lines[0]
    natoms = int(lines[1])
    atoms = lines[2:2+natoms]
    box = list(map(float, lines[2+natoms].split()))[:3]  # ✅ 修复：只保留 X Y Z
    return title, natoms, atoms, np.array(box)


def compute_center(atoms):
    coords = []
    for line in atoms:
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        coords.append([x, y, z])
    return np.mean(coords, axis=0)

def shift_atoms(atoms, shift_vec):
    new_atoms = []
    for line in atoms:
        x = float(line[20:28]) + shift_vec[0]
        y = float(line[28:36]) + shift_vec[1]
        z = float(line[36:44]) + shift_vec[2]
        newline = line[:20] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[44:]
        new_atoms.append(newline)
    return new_atoms

def write_gro(filename, title, natoms, atoms, box):
    with open(filename, 'w') as f:
        f.write(title)
        f.write(str(natoms) + '\n')
        f.writelines(atoms)
        f.write(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("❌ 用法错误：请提供一个 .gro 文件名作为输入。\n示例：python center_gro_precisely.py complex.gro")
        sys.exit(1)

    input_file = sys.argv[1]
    if not os.path.isfile(input_file):
        print(f"❌ 输入文件不存在：{input_file}")
        sys.exit(1)

    output_file = os.path.splitext(input_file)[0] + "_centered.gro"

    title, natoms, atoms, box = read_gro(input_file)
    center_old = compute_center(atoms)
    center_new = box / 2
    shift = center_new - center_old
    new_atoms = shift_atoms(atoms, shift)
    write_gro(output_file, title, natoms, new_atoms, box)

    print(f"✅ 已将复合体从 {center_old.round(3)} 平移到盒子中心 {center_new.round(3)}")
    print(f"👉 新文件输出为：{output_file}")
