#!/usr/bin/env python3
"""trace_downstream.py — 扫描给定 h5ad 版本的下游依赖链

用法:
    python scripts/trace_downstream.py results/04_embedded_v1.h5ad

输出每个直接/间接消费该文件的 notebook + 其 PARAMS 中的 UPSTREAM_PATH。
帮助 PI 判断：改了某个 stage 后，哪些下游需要重跑。
"""
import sys
import os
import json
import re

def find_notebooks(root="notebooks"):
    """递归找所有 .ipynb"""
    nbs = []
    for dirpath, _, filenames in os.walk(root):
        for f in filenames:
            if f.endswith(".ipynb") and not f.startswith("."):
                nbs.append(os.path.join(dirpath, f))
    return sorted(nbs)

def extract_upstream(nb_path):
    """从 notebook 的 PARAMS cell 中提取 UPSTREAM_PATH"""
    with open(nb_path, "r") as f:
        nb = json.load(f)

    upstreams = []
    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        source = "".join(cell.get("source", []))
        # 匹配 UPSTREAM_PATH = "..."
        matches = re.findall(r'UPSTREAM_PATH\s*=\s*["\']([^"\']+)["\']', source)
        upstreams.extend(matches)
    return upstreams

def trace(target_file, notebooks):
    """找所有直接或间接依赖 target_file 的 notebook"""
    # 第一层：直接消费
    direct = []
    dep_map = {}  # nb_path -> [upstream_paths]

    for nb in notebooks:
        upstreams = extract_upstream(nb)
        dep_map[nb] = upstreams
        if target_file in upstreams or os.path.basename(target_file) in [os.path.basename(u) for u in upstreams]:
            direct.append(nb)

    # 递归：间接消费（如 04 → 05 → 06，改 04 后 06 也需重跑）
    all_affected = set(direct)
    queue = list(direct)
    while queue:
        current = queue.pop(0)
        # 找 current notebook 的 OUTPUT_PATH
        with open(current, "r") as f:
            nb = json.load(f)
        for cell in nb.get("cells", []):
            if cell.get("cell_type") != "code":
                continue
            source = "".join(cell.get("source", []))
            outputs = re.findall(r'OUTPUT_PATH\s*=\s*["\']([^"\']+)["\']', source)
            for out in outputs:
                # 找谁消费了这个 output
                for nb2 in notebooks:
                    if nb2 in all_affected:
                        continue
                    if out in dep_map.get(nb2, []) or os.path.basename(out) in [os.path.basename(u) for u in dep_map.get(nb2, [])]:
                        all_affected.add(nb2)
                        queue.append(nb2)

    return direct, all_affected - set(direct)

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    target = sys.argv[1]
    notebooks = find_notebooks()

    print(f"目标文件: {target}")
    print(f"扫描 {len(notebooks)} 个 notebook...\n")

    direct, indirect = trace(target, notebooks)

    if direct:
        print("直接消费（UPSTREAM_PATH 指向此文件）:")
        for nb in direct:
            print(f"  → {nb}")
    else:
        print("无直接消费者")

    if indirect:
        print(f"\n间接消费（通过依赖链传递）:")
        for nb in indirect:
            print(f"  ⤳ {nb}")

    if direct or indirect:
        print(f"\n总计需重跑: {len(direct) + len(indirect)} 个 notebook")
        print("建议按编号从小到大依次重跑。")
    else:
        print("\n该文件无下游消费者——可安全修改不影响其他 stage。")

if __name__ == "__main__":
    main()
