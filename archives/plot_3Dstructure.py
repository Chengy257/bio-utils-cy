import os
import glob
import argparse
from pymol import cmd

HYDROPHOBICITY = {
    'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5,
    'CYS': 2.5, 'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4,
    'HIS': -3.2, 'ILE': 4.5, 'LEU': 3.8, 'LYS': -3.9,
    'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6, 'SER': -0.8,
    'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
}

def setup_hydrophobicity_colors():
    """定义更柔和、论文友好的疏水性配色"""
    cmd.set_color("hydro_blue",  [0.30, 0.55, 0.85])  # 柔和蓝
    cmd.set_color("hydro_white", [1.00, 1.00, 1.00])  # 中性
    cmd.set_color("hydro_orange",[0.95, 0.55, 0.20])  # 柔和橙

def initialize_pymol():
    """初始化PyMOL以无GUI模式运行"""
    cmd.reinitialize()

def color_by_secondary_structure(obj_name):
    """按二级结构着色（helix / sheet / loop）"""

    # 计算二级结构
    cmd.dss(obj_name)

    cmd.show_as("cartoon", obj_name)

    # 定义论文友好颜色
    cmd.set_color("ss_helix", [0.20, 0.35, 0.70])   # 深蓝
    cmd.set_color("ss_sheet", [0.90, 0.60, 0.20])  # 金橙
    cmd.set_color("ss_loop",  [0.75, 0.75, 0.75])  # 浅灰

    # 按结构类型着色
    cmd.color("ss_helix", f"{obj_name} and ss H")
    cmd.color("ss_sheet", f"{obj_name} and ss S")
    cmd.color("ss_loop",  f"{obj_name} and ss L")
    cmd.color("ss_loop",  f"{obj_name} and ss ''")

    # 不再使用 putty（避免新版本歧义）
    cmd.set("cartoon_transparency", 0.0)


def color_by_hydrophobicity(obj_name):
    """根据疏水性值着色（优化配色）"""

    # 写入疏水性到 b-factor
    for resn, value in HYDROPHOBICITY.items():
        cmd.alter(f"{obj_name} and resn {resn}", f"b = {value}")

    # 固定范围，避免被极端值拉爆
    cmd.spectrum(
        "b",
        "hydro_blue hydro_white hydro_orange",
        obj_name,
        minimum=-4.5,
        maximum=4.5
    )

    cmd.show_as("cartoon", obj_name)

    # cartoon 微调，更“干净”
    cmd.set("cartoon_putty", 1)
    cmd.set("cartoon_putty_scale_min", 0.6)
    cmd.set("cartoon_putty_scale_max", 1.4)
    cmd.set("cartoon_putty_radius", 0.4)


def generate_images(input_folder, output_folder, width, height, dpi):
    """处理所有结构文件生成PNG图片"""
    for struct_file in glob.glob(os.path.join(input_folder, "**/*.*"), recursive=True):
        if not struct_file.lower().endswith(('.pdb', '.cif')):
            continue
        
        cmd.delete("all")
        cmd.load(struct_file)
        obj_name = cmd.get_object_list()[0]
        # color_by_hydrophobicity(obj_name)        
        color_by_secondary_structure(obj_name)
        cmd.zoom()
        cmd.orient()
        
        rel_path = os.path.relpath(struct_file, input_folder)
        base_name = os.path.splitext(rel_path)[0]
        png_path = os.path.join(output_folder, f"{base_name}.png")
        os.makedirs(os.path.dirname(png_path), exist_ok=True)
        
        cmd.png(png_path, width=width, height=height, dpi=dpi, ray=1)
        print(f"已生成：{png_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="批量生成结构文件的PNG图片（疏水性着色）",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", help="输入文件夹路径（包含PDB/CIF文件）")
    parser.add_argument("-o", "--output", default="output", help="输出文件夹路径")
    parser.add_argument("--width", type=int, default=2400, help="图片宽度（像素）")
    parser.add_argument("--height", type=int, default=2400, help="图片高度（像素）")
    parser.add_argument("--dpi", type=int, default=600, help="输出分辨率")
    args = parser.parse_args()

    initialize_pymol()
    os.makedirs(args.output, exist_ok=True)
    generate_images(args.input, args.output, args.width, args.height, args.dpi)
    print(f"\n处理完成！结果已保存至：{os.path.abspath(args.output)}")
