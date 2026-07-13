---
title: "Manifest Schema 迁移顺序"
tags:
  - scRNA-seq
  - manifest
  - schema
  - migration
created: "2026-07-13"
updated: "2026-07-13"
---

# Manifest Schema 迁移顺序

## 目标与边界

本迁移把当前运行 manifest 逐步升级为可严格验证的 v1 schema，同时保持每个中间 PR 可审查、可回退。迁移只处理运行契约、provenance 和 artifact 语义；PR2 的 counts 来源、SoupX 使用范围和数据尺度语义不得混入本序列。

## 实施顺序

### g1a：仅定义 v1 schema

- 独立加入 v1 schema、字段约束和纯验证测试。
- 不改任何 notebook、runner、resume、promotion 或 checkpoint 的运行行为。
- 不启用严格校验，不让未迁移 producer 提前失败。

### g2：补齐共享 provenance 与 artifact 能力

- 时间统一为 UTC，避免本地时区产生不可比较记录。
- 增加受限文件树 fingerprint，用于记录输入目录的确定性状态。
- 记录 R 运行环境 provenance。
- 提供统一 artifact 构造 helper，集中生成受管路径、角色、大小和哈希字段。

### g3-g14：逐 producer 迁移

按小 PR 依次迁移四个 Stage 01 notebook、Stage 01 template，以及 Stage 02、03、04、05、06、06b、06c。每个 producer 必须：

- 写出符合 v1 schema 的 manifest，并登记全部 artifacts，而非只登记 primary output。
- 记录实际使用的输入和 fingerprint，不以配置值代替运行期真实输入。
- 保留现有科研参数与 notebook 可见性，不借 schema 迁移改变分析方法。
- 让 06b 完整登记报告 payload，包括正文、图片和必要辅助文件。

迁移期间，`allow_legacy=True` 只允许用于读取旧 manifest 的过渡兼容；它不得生成、改写或晋升 legacy manifest。

### g15：统一激活严格运行校验

仅在 g3-g14 的 producer inventory 全部迁移后，统一启用严格的 `validate_checkpoint`、artifact 校验、resume、promotion 和 runner 判定。激活时必须：

- checkpoint 与全部 artifacts 同时通过路径、角色和哈希验证。
- 将已验证 manifest 的哈希绑定到后续 promotion，在最终移动前复核，阻断 manifest TOCTOU 替换。
- 从本阶段起，promotion 绝不得使用 `allow_legacy=True`；legacy manifest 只能读取和诊断，不能晋升。

### g16：退役活动 legacy fallback

完成 producer inventory，确认所有活动 producer 的 legacy 输出为零后，删除运行路径中的 legacy fallback。历史 manifest 的离线读取可保留为明确隔离的兼容入口，但不得参与 resume、promotion、runner 成功判定或新 checkpoint 生成。

## 合并门槛

每个 PR 独立通过 schema/行为测试和静态检查；g15 前不得提前激活严格模式，g16 前必须提供活动 producer inventory 为零的证据。
