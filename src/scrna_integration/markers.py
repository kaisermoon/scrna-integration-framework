"""标记物库加载器：读取 ``references/markers/{组织}_{用途}.csv``。

ADR-0005 论证：load+filter+groupby 这套操作在 stage 6 标注、
逐簇深度分析和基因集评分 notebook 中反复出现。集中处理 *role* 语义
（canonical/optional/negative）防止误用。

面向非计算机专业 PI/学生：只需一个调用即可获取按细胞类型组织的标记物列表。
"""

from __future__ import annotations

import pandas as pd


def load_markers(
    csv_path: str,
    roles: tuple[str, ...] | None = ("canonical", "optional"),
) -> dict[str, list[str]] | dict[str, dict[str, list[str]]]:
    """加载标记物 CSV 并返回 ``细胞类型 → 标记物列表``（按 role 过滤）。

    面向非 CS 学生：CSV 是课题组长期积累的标记物知识库，每行一个基因。
    本函数帮你按 cell_type 分组并按 role 筛选——避免每次打开 notebook
    都要重写 pandas 的 groupby+filter 操作。

    Args:
        csv_path: 任意文件系统路径（相对或绝对）。
        roles: role 过滤元组。默认 ``("canonical", "optional")``
            返回基因集评分 / dotplot 的常用情形。
            ``roles=("negative",)`` 仅返回阴性标记物。
            ``roles=None`` 返回完整三层字典
            ``{细胞类型: {canonical: [...], optional: [...], negative: [...]}}``。

    Returns:
        以细胞类型为 key 的字典——*roles* 为元组时是扁平列表，
        *roles* 为 ``None`` 时是嵌套字典。

    CSV 格式: ``tissue, cell_type, marker, role, reference, notes``

    典型用法:
        markers = load_markers("references/markers/gastric_epithelial.csv")
        # → {"SPEM": ["TFF2", "MUC6", ...], "pit_cell": ["MUC5AC", ...]}
        sc.pl.dotplot(adata, var_names=markers, groupby="leiden")
    """
    df = pd.read_csv(csv_path)

    if roles is None:
        # 完整三层模式：{cell_type: {canonical: [...], optional: [...], negative: [...]}}
        result: dict[str, dict[str, list[str]]] = {}
        for role in ("canonical", "optional", "negative"):
            role_df = df[df["role"] == role]
            for ct, group in role_df.groupby("cell_type"):
                result.setdefault(ct, {}).setdefault(role, [])
                result[ct][role] = group["marker"].tolist()
        # 确保每个 cell_type 都有三种 role 键（没有的为空列表）
        for ct in result:
            for role in ("canonical", "optional", "negative"):
                result[ct].setdefault(role, [])
        return result

    # 扁平模式：按 role 过滤，返回 {cell_type: [markers]}
    filtered = df[df["role"].isin(roles)]
    return {
        ct: group["marker"].tolist() for ct, group in filtered.groupby("cell_type")
    }
