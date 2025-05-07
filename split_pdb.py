
import sys

def split_pdb_by_TER(input_pdb_path):
    protein_lines = []
    ligand_lines = []
    connect_lines = []

    current_block = []
    is_protein = True  # 初始设定第一个为蛋白质

    with open(input_pdb_path, 'r') as file:
        for line in file:
            if line.startswith("CONECT"):
                connect_lines.append(line)
                continue

            if line.startswith("TER"):
                if is_protein:
                    protein_lines += current_block + [line]
                    is_protein = False
                else:
                    ligand_lines += current_block + [line]
                current_block = []
            elif line.startswith(("ATOM", "HETATM", "ANISOU", "SIGUIJ", "SIGATM")):
                current_block.append(line)
            else:
                continue

    if current_block:
        if is_protein:
            protein_lines += current_block
        else:
            ligand_lines += current_block

    ligand_lines += connect_lines

    base = input_pdb_path.rsplit(".", 1)[0]
    protein_path = "protein.pdb"
    ligand_path = "ligand.pdb"

    with open(protein_path, 'w') as f:
        f.writelines(protein_lines)
    with open(ligand_path, 'w') as f:
        f.writelines(ligand_lines)

    print(f"✅ 拆分完成：\n蛋白质部分：{protein_path}\n小分子部分：{ligand_path}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("用法：python split_pdb.py 输入文件.pdb")
        sys.exit(1)
    split_pdb_by_TER(sys.argv[1])
