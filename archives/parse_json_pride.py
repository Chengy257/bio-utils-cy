#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import json
import ast
import argparse
from pathlib import Path
import glob
import pandas as pd

def load_json(file_path):
    """兼容 PRIDE JSON 格式，返回 Python dict"""
    text = Path(file_path).read_text(encoding="utf-8")
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        # 有些导出可能是 Python dict 字符串而非严格 JSON
        return ast.literal_eval(text)

def _first_nonempty(*vals):
    for v in vals:
        if v is None:
            continue
        if isinstance(v, str) and v.strip() == "":
            continue
        return v
    return None

def extract_instruments(data):
    """提取 instruments 的易读格式（优先 value -> name -> accession -> cvLabel），并去重，返回逗号分隔字符串"""
    insts = data.get("instruments", []) or []
    models = []
    for inst in insts:
        # 优先取 value（通常是真正的型号），fallback 为 name / accession / cvLabel / @type
        model = _first_nonempty(inst.get("value"), inst.get("name"), inst.get("accession"), inst.get("cvLabel"), inst.get("@type"))
        if model:
            model = str(model).strip()
            if model not in models:
                models.append(model)
    return ", ".join(models)

def extract_core_info(data):
    """提取核心元信息（包含改进后的 instruments）"""
    return {
        "Accession": data.get("accession"),
        "Title": data.get("title"),
        "Description": data.get("projectDescription"),
        "DOI": data.get("doi"),
        "Submission Date": data.get("submissionDate"),
        "Publication Date": data.get("publicationDate"),
        "License": data.get("license"),
        "Organisms": ", ".join([org.get("name","") for org in data.get("organisms",[])]),
        "Organism Parts": ", ".join([p.get("name","") for p in data.get("organismParts",[])]),
        # 使用改进函数
        "Instruments": extract_instruments(data),
        "Softwares": ", ".join([s.get("name","") for s in data.get("softwares",[])]),
        "Experiment Types": ", ".join([e.get("name","") for e in data.get("experimentTypes",[])]),
        "Submitters": ", ".join([s.get("name","") for s in data.get("submitters",[])]),
        "Affiliations": ", ".join(data.get("affiliations",[])),
        "Countries": ", ".join(data.get("countries",[])),
        "Keywords": ", ".join(data.get("keywords",[])),
        "References": "; ".join([r.get("referenceLine","") for r in data.get("references",[])]),
        "Identified PTMs": ", ".join([ptm.get("name","") for ptm in data.get("identifiedPTMStrings",[])]),
    }

def gather_files(input_paths, recursive=True):
    files = []
    for p in input_paths:
        # expand user and handle glob patterns
        if any(ch in p for ch in ["*", "?", "["]):
            expanded = glob.glob(p, recursive=recursive)
            files.extend([f for f in expanded if Path(f).is_file()])
            continue

        path = Path(p).expanduser()
        if path.is_dir():
            if recursive:
                files.extend([str(f) for f in path.rglob("*.json")])
            else:
                files.extend([str(f) for f in path.glob("*.json")])
        elif path.is_file():
            files.append(str(path))
        else:
            print(f"[WARN] 未找到路径或模式：{p}")
    # 去重且按路径排序
    files = sorted(list(dict.fromkeys(files)))
    return files

def main(args):
    files = gather_files(args.inputs, recursive=not args.no_recursive)
    if not files:
        print("没有找到任何 JSON 文件。请检查路径/模式。")
        return

    records = []
    errors = []
    for f in files:
        try:
            data = load_json(f)
            rec = extract_core_info(data)
            rec["Source File"] = Path(f).name
            records.append(rec)
        except Exception as e:
            errors.append((f, str(e)))
            print(f"[ERROR] 解析失败: {f} -> {e}")

    if not records:
        print("没有成功解析的记录。")
        return

    df = pd.DataFrame(records)

    # 输出 Excel 和/或 CSV
    out_base = Path(args.output)
    suffix = out_base.suffix.lower()
    if args.format in ("xlsx","both"):
        xlsx_path = out_base if suffix in (".xlsx",) else out_base.with_suffix(".xlsx")
        df.to_excel(xlsx_path, index=False, engine="openpyxl")
        print(f"✅ 已输出 Excel: {xlsx_path}")
    if args.format in ("csv","both"):
        csv_path = out_base if suffix in (".csv",) else out_base.with_suffix(".csv")
        # 使用 utf-8-sig 以便在 Excel 中正确识别中文
        df.to_csv(csv_path, index=False, encoding="utf-8-sig")
        print(f"✅ 已输出 CSV: {csv_path}")

    print(f"共处理文件: {len(files)}，成功: {len(records)}，失败: {len(errors)}")
    if errors:
        print("失败列表 (前 10 条):")
        for f, e in errors[:10]:
            print(" -", f, ":", e)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="批量解析 PRIDE JSON 元数据，优先提取 instruments.name/value")
    parser.add_argument("inputs", nargs="+", help="文件/目录/glob 模式 (例如 data/*.json 或 ./folder)")
    parser.add_argument("-o", "--output", default="pride_metadata.xlsx", help="输出文件名或基名 (默认: pride_metadata.xlsx)")
    parser.add_argument("-f", "--format", choices=["xlsx","csv","both"], default="xlsx", help="输出格式 (默认 xlsx)")
    parser.add_argument("--no-recursive", action="store_true", help="遍历目录时不递归（只扫描顶层）")
    args = parser.parse_args()
    main(args)
