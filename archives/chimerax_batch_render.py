#!/usr/bin/env python3
# chimerax_batch_render.py

import os
import glob
import subprocess
import tempfile

def create_chimerax_script(pdb_files, output_folder):
    """创建ChimeraX脚本内容"""
    script_content = """
# 设置背景和光照
set bgColor white
lighting soft
graphics silhouettes true

runscript inline
import os

output_folder = "{}"
pdb_files = {}

for pdb_file in pdb_files:
    pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
    print(f"Processing: {{pdb_name}}")
    
    try:
        # 加载PDB文件
        run(f"open {{pdb_file}}")
        
        # 设置显示样式
        run("hide atoms")
        run("show cartoons")
        run("color bychain")
        
        # 优化视图
        run("view")
        run("lighting soft")
        
        # 保存图片
        output_path = os.path.join(output_folder, f"{{pdb_name}}.png")
        run(f"save {{output_path}} width 600 height 600 supersample 2")
        
        # 关闭模型
        run("close")
        print(f"Saved: {{output_path}}")
        
    except Exception as e:
        print(f"Error processing {{pdb_file}}: {{e}}")
        run("close")

print("批量渲染完成！")
endscript
""".format(output_folder, pdb_files)
    
    return script_content

def batch_render_with_chimerax(pdb_folder, output_folder):
    """使用ChimeraX批量渲染PDB文件"""
    
    # 创建输出文件夹
    os.makedirs(output_folder, exist_ok=True)
    
    # 获取所有PDB文件
    pdb_files = glob.glob(os.path.join(pdb_folder, "*.pdb"))
    
    if not pdb_files:
        print(f"在 {pdb_folder} 中未找到PDB文件")
        return
    
    print(f"找到 {len(pdb_files)} 个PDB文件")
    
    # 创建临时脚本文件
    script_content = create_chimerax_script(pdb_files, output_folder)
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.cxc', delete=False) as f:
        f.write(script_content)
        script_path = f.name
    
    try:
        # 运行ChimeraX
        cmd = ['chimerax', '--nogui', '--exit', script_path]
        print("启动ChimeraX进行批量渲染...")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            print("ChimeraX渲染完成！")
        else:
            print(f"ChimeraX运行出错: {result.stderr}")
            
    finally:
        # 清理临时文件
        os.unlink(script_path)

if __name__ == "__main__":
    pdb_folder = "./pdb_files"
    output_folder = "./structure_images"
    batch_render_with_chimerax(pdb_folder, output_folder)
