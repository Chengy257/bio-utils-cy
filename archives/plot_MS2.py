import argparse
import matplotlib.pyplot as plt
from lxml import etree
import pymzml
import re
import numpy as np
from pyteomics import mass, mgf

# 配置参数解析
parser = argparse.ArgumentParser(description='MSGF+结果可视化工具')
parser.add_argument('--mzid', type=str, required=True, help='MSGF+的mzid文件路径')
parser.add_argument('--spectra', type=str, required=True, help='质谱文件路径（支持.mzML或.mgf）')
parser.add_argument('--protein', type=str, required=True, help='目标蛋白ID')
parser.add_argument('--output', type=str, default='output.png', help='输出图片路径')
args = parser.parse_args()

# 常量定义
NS = {'ns': 'http://psidev.info/psi/pi/mzIdentML/1.2'}
ION_OFFSET = 0.1  # 标注偏移量
COLORS = {'b': '#FF3030', 'y': '#4876FF'}

def get_namespaces(tree):
    namespaces = {}
    for elem in tree.iter():
        if elem.tag.startswith('{'):
            ns = elem.tag.split('}')[0][1:]
            namespaces['ns'] = ns
            break
    return namespaces

# 解析mzid文件（兼容蛋白ID模糊匹配）
def parse_mzid(target_protein_id):
    psm_list = []
    tree = etree.parse(args.mzid)
    NS = get_namespaces(tree)
    
    # 清洗目标蛋白ID
    clean_target = re.split(r'[\.\s_]', target_protein_id)[0].upper()
    print(f"[DEBUG] 清洗后的目标蛋白ID: {clean_target}")
    
    for spectrum in tree.xpath('//ns:SpectrumIdentificationResult', namespaces=NS):
        spectrum_id = spectrum.get('spectrumID')
        scan_num = int(re.search(r'scan=(\d+)', spectrum_id).group(1)) if 'scan=' in spectrum_id else int(spectrum_id)
        
        for item in spectrum.xpath('.//ns:SpectrumIdentificationItem', namespaces=NS):
            # 提取所有可能的蛋白标识符
            protein_refs = item.xpath('''
                .//ns:ProteinAccession/text() |
                .//ns:DBSequence/@accession |
                .//ns:ProteinDetectionHypothesis/@id
            ''', namespaces=NS)
            
            # 模糊匹配逻辑
            matched = False
            for ref in protein_refs:
                clean_ref = re.split(r'[\.\s_]', ref)[0].upper()
                if clean_ref == clean_target:
                    matched = True
                    break
            
            if not matched:
                continue
                
            # 提取肽段信息...
            psm_list.append({...})
    
    print(f"[DEBUG] 找到的PSM列表: {psm_list}")
    return psm_list

# 解析带修饰的肽段序列（同上）
def parse_peptide_sequence(peptide_ref):
    base_seq = re.sub(r'\[.*?\]', '', peptide_ref.split('#')[0])
    return base_seq.upper()

# 计算理论碎片离子（同上）
def calculate_fragments(sequence, charge, modifications):
    calc = mass.Composition(sequence=sequence)
    frags = {'b': [], 'y': []}
    mod_dict = {pos: delta for pos, delta in modifications}
    
    for i in range(1, len(sequence)):
        # 计算b离子
        b_comp = calc[:i]
        for pos in mod_dict:
            if pos <= i:
                b_comp['mass'] += mod_dict[pos]
        frags['b'].append(mass.calculate_mass(b_comp, charge=1))
        
        # 计算y离子
        y_comp = calc[-i:]
        for pos in mod_dict:
            if pos > (len(sequence) - i):
                y_comp['mass'] += mod_dict[pos]
        frags['y'].append(mass.calculate_mass(y_comp, charge=1))
    
    return frags

# 加载质谱数据（新增MGF支持）
def load_spectra(file_path):
    if file_path.lower().endswith('.mgf'):
        return {spec['params']['scans']: spec for spec in mgf.read(file_path)}
    elif file_path.lower().endswith(('.mzml', '.mzxml')):
        return pymzml.run.Reader(file_path)
    else:
        raise ValueError("不支持的质谱文件格式")

# 主流程
if __name__ == "__main__":
    # 步骤1: 解析mzid获取目标PSM
    psms = parse_mzid(args.protein)
    if not psms:
        raise ValueError(f"未找到蛋白 {args.protein} 的PSM信息")

    # 步骤2: 加载质谱数据（兼容MGF/mzML）
    spectra_data = load_spectra(args.spectra)

    # 步骤3: 遍历每个PSM绘制谱图
    for idx, psm in enumerate(psms):
        # 获取谱图数据
        if isinstance(spectra_data, dict):  # MGF模式
            spec = spectra_data.get(str(psm['scan']))
            if not spec:
                print(f"[WARNING] 在MGF文件中未找到scan={psm['scan']}")
                continue
            mz = spec['m/z array']
            intensity = spec['intensity array']
        else:  # mzML模式
            spec = spectra_data[str(psm['scan'])]
            if spec.ms_level != 2:
                continue
            mz = spec.mz
            intensity = spec.i
        
        # 计算理论碎片
        frags = calculate_fragments(psm['sequence'], psm['charge'], psm['mods'])
        
        # 创建画布
        plt.figure(figsize=(12, 6))
        plt.title(f"Scan {psm['scan']} | {psm['sequence']}\nProtein: {args.protein}", fontsize=12)
        
        # 绘制实际谱图
        plt.stem(mz, intensity, linefmt='gray', markerfmt=' ', basefmt=' ')
        
        # 标注碎片离子
        for ion_type in ['b', 'y']:
            for i, theo_mz in enumerate(frags[ion_type]):
                nearest_idx = np.argmin(np.abs(mz - theo_mz))
                if abs(mz[nearest_idx] - theo_mz) < 0.1:
                    label = f"{ion_type}{i+1}"
                    plt.vlines(theo_mz, 0, intensity[nearest_idx], 
                              colors=COLORS[ion_type], 
                              linestyles='dashed',
                              alpha=0.7)
                    plt.text(theo_mz + ION_OFFSET, 
                            intensity[nearest_idx], 
                            label,
                            color=COLORS[ion_type],
                            fontsize=8)
        
        # 显示氨基酸序列
        seq_text = ' '.join(list(psm['sequence']))
        plt.figtext(0.5, 0.01, seq_text, 
                   ha='center', 
                   fontsize=10,
                   bbox=dict(facecolor='lightgray', alpha=0.5))
        
        # 保存结果
        output_path = f"{args.output}_{idx+1}.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()