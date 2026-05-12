#!/usr/bin/env python3
# batch_alphafold_pymol.py

import os
import glob
import sys
import time

sys.path.append('/usr/lib/python3/dist-packages/')

import __main__
__main__.pymol_argv = ['pymol', '-qc']

import pymol
from pymol import cmd

def setup_alphafold_colors():
    """设置AlphaFold官方颜色"""
    # 定义pLDDT颜色
    cmd.set_color('af_very_high', [0.051, 0.341, 0.827])  # 深蓝色 >90
    cmd.set_color('af_confident', [0.416, 0.796, 0.945])  # 浅蓝色 70-90
    cmd.set_color('af_low', [0.996, 0.851, 0.212])       # 黄色 50-70
    cmd.set_color('af_very_low', [0.992, 0.490, 0.302])  # 橙色 <50

def color_by_plddt(pdb_name):
    """根据pLDDT值着色"""
    try:
        # 应用AlphaFold颜色方案
        cmd.color('af_very_high', f'{pdb_name} and b > 90')
        cmd.color('af_confident', f'{pdb_name} and b > 70 and b <= 90')
        cmd.color('af_low', f'{pdb_name} and b > 50 and b <= 70') 
        cmd.color('af_very_low', f'{pdb_name} and b <= 50')
    except:
        # 备选方案：使用spectrum
        # cmd.spectrum('b', 'tv_red_yellow_green_cyan_marine', pdb_name, minimum=50, maximum=90)]
        None

def batch_generate_alphafold_images(pdb_folder, output_folder):
    """批量生成AlphaFold结构图片，使用pLDDT着色"""
    
    os.makedirs(output_folder, exist_ok=True)
    pdb_files = glob.glob(os.path.join(pdb_folder, "*.pdb"))
    
    if not pdb_files:
        print(f"未在 {pdb_folder} 中找到PDB文件")
        return
    
    pymol.finish_launching()
    time.sleep(1)
    
    # 设置AlphaFold颜色
    setup_alphafold_colors()
    
    for pdb_file in pdb_files:
        try:
            pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
            print(f"处理: {pdb_name}")
            
            cmd.delete('all')
            cmd.load(pdb_file, pdb_name)
            time.sleep(0.5)
            
            # 设置显示样式
            cmd.hide('everything')
            cmd.show('cartoon', pdb_name)
            
            # 应用pLDDT着色
            color_by_plddt(pdb_name)
            
            # 设置背景和渲染
            cmd.bg_color('white')
            cmd.set('antialias', 1)
            cmd.center(pdb_name)
            cmd.orient(pdb_name)
            cmd.zoom(pdb_name, buffer=2)
            cmd.refresh()
            time.sleep(0.2)
            
            # 生成图片
            output_path = os.path.join(output_folder, f"{pdb_name}_plddt.png")
            cmd.png(output_path, width=800, height=600, dpi=150, ray=1)
            
            print(f"已生成: {output_path}")
            
        except Exception as e:
            print(f"处理 {pdb_file} 时出错: {str(e)}")
            continue
    
    cmd.quit()
    print("AlphaFold批量着色完成！")

if __name__ == "__main__":
    pdb_folder = "./pdb_files"
    output_folder = "./alphafold_images"
    batch_generate_alphafold_images(pdb_folder, output_folder)
