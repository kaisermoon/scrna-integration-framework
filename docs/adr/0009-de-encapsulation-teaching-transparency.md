---
status: accepted
---

# 从薄框架进一步回退到「教学透明」：删 sweep、摊平 IO 琐碎 helper、注释中文化

第三轮 grilling 后框架已是"3 函数 + 约定"的薄壳，但 PI 在进入实施期、审阅已合并的 PR-0a~3b 代码后给出更强的方向纠偏：**对非计算机专业的 PI 与学生，当前封装仍过度，学习成本过高**。本 ADR 记录由此触发的进一步回退，以及它与 ADR-0001/0003/0004（薄框架/朴素优先/删除清单）的关系——它不是推翻那些 ADR，而是把"薄"的判据从"框架作者觉得薄"收紧为"非 CS 学生打开 notebook 能逐行看懂"。

## 触发的事实（Explore 量化）

- `io.py` 734 行 / 24 个私有 helper，其中约 250 行是"简单操作被包一层"：`_propagate_disease_system`（1 行赋值）、`_enforce_species`（2 行）、`_read_h5ad`/`_read_h5`（各 3 行）、`_record_raw_path`、`_generate_cell_id`、`_apply_obs_mapping`、`_inject_constants`、`_rename_original_annotations` 等。
- `read_with_manifest` 在 notebook 里是 1 行黑盒，学生看不到 cell_id 生成、QC 计算、基因 ID 同步等底层流程。
- `sweep(fn, adata, candidates, scorer, output_dir)` 把一个简单 for 循环抽象成 5 参数 + 回调（`fn=`/`scorer=`）。非 CS 学生惧怕回调，却能轻松理解 for 循环 + `adata.copy()` + `sc.pp.*` + `pd.DataFrame`。
- 全部注释为英文，且多解释 what 不解释 why——对已懂流程的人够用，对初学者学习成本高。

这与项目原始构思第 61-63 行的明确要求一致："面向 PI 和非计算机专业的学生，代码一定要简洁易懂，有充分的注释……尽可能使用 scanpy 原生函数，不要进行过分的函数封装或新建一套类的体系。"

## 决定

1. **删除 `sweep` 框架函数**（及其 `_write_report`）。stage4/stage5 notebook 改为**显式 for 循环**，循环体内直接调 scanpy 原生 API + 直接调 `scorers` 里的指标函数。`sweep` 退出框架表面。
2. **`scorers` 保留为"可直接调用的指标函数"**（不是回调）。notebook 在 for 循环里 `m = integration_metrics(adata_copy)` 这样直接调用——透明、无 `fn=`/`scorer=` 抽象。scIB 套件（kBET/iLISI 等）学生不会自己重写，留作可直接调用的 helper 有真实价值。
3. **摊平 `io.py` 的琐碎 helper**：约 10 个 1-3 行的 wrapper 内联回 `read_with_manifest` 主体，使其从上到下像带中文注释的教程读下来。**保留**真正复杂的（多源 mtx 目录发现、基因 ID 双向同步 mygene、临床表 join、manifest 校验）——内联进 notebook 反而更乱，但这些保留的 helper 也要补充中文注释讲清在做什么。
4. **注释全面中文化且充分**：所有 notebook + src 的注释改中文，关键步骤解释 why（为什么 target_sum=1e4、为什么 seurat_v3 HVG flavor、为什么这样同步基因 ID），不只 what。
5. 框架表面从"3 函数"收缩为 **`read_with_manifest` + `load_markers` + 可直接调用的 `scorers` 指标函数**。CONTEXT 的 Sweep 词条删除。

## Considered Options

- **保留 sweep 仅加中文注释**：拒绝。回调抽象本身是学习成本的来源，注释救不了"学生看到 `fn=`/`scorer=` 就却步"。PI 明确选删。
- **连 scorers 一起删、全内联进 notebook**：拒绝（暂）。scIB 指标套件是真实复杂度，每个 notebook 内联会产生大段学生看不懂也不会维护的代码；保留为"可直接调用、无回调"的函数是封装与透明的平衡点。若后续仍嫌重可再议。
- **把 read_with_manifest 整个拆进 notebook**：拒绝。多源异构 IO 是框架存在的核心理由（原始构思第一诉求），拆进 notebook 会让每个分析重复几百行脆弱的格式处理。保留为单函数，但主体线性化 + 中文注释 + notebook 配"它做了什么"说明 cell。

## Consequences

- PR-2 引入的 `sweep` 被删；stage4/stage5 notebook（PR-3a/3b 已合并）改写为显式循环。这是有意推翻 PR-2 的一部分，记录在案。
- `scrna_integration` 的 `__init__.py` re-export 改为 `read_with_manifest` + `load_markers`（+ 可选 `scorers` 模块导入）。
- SPEC 的 "The Three Functions" 节改为 "Two Functions + Scorers"；所有 sweep 调用示例改为显式循环；stage4/5 cell 序列更新。
- code-reviewer 今后的"过度封装"判据收紧：新增任何 helper/抽象，问"非 CS 学生打开 notebook 能否逐行看懂"，而非只问"是否填了 scanpy 空白"。
- 注释语言从英文改中文是项目级约定（科研代码面向中文 PI + 学生），写入 SPEC。
