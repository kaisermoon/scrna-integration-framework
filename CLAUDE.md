---
title: "scRNA-seq整合分析框架 项目级指令"
type: project-instructions
project_id: "scrna-integration-framework"
created: "2026-06-08"
updated: "2026-06-10"
---

# CLAUDE.md — scRNA-seq整合分析框架

> 本文件是**项目级**指令，叠加在顶级 `~/AI-OS/CLAUDE.md` + `SOUL.md` 之上。进入本项目时自动加载（CWD 触发）。
> 行为约束指针还包括：`项目/_GitHub项目规范.md`（GitHub 仓库项目规则，进入时主 Agent 主动 Read）。

## 一、最高纪律：默认委派 subagent，主 Agent 只做思考/决策/调度，只碰"摘要 + 决策点 + 路径"

**本项目反复出现主 Agent 越界自己干执行活、以及把长原文整篇读进来的问题。强制纠偏：判据不是"这是思考还是执行"，而是"这个动作会不会把大量原文 / 数据 / 日志拉进主 Agent 窗口"。会的，一律委派或定向化。**

主 Agent **只做**：决策与判断、调度、写 brief、读 subagent 回报、verdict 路由、size/红线裁决、终审、ADR 撰写、与 PI 对话。

主 Agent **不亲自做**（必须委派）：写/改代码与 notebook、跑 nbconvert/pytest/ruff、`gh pr merge`/worktree 清理/分支删除、CI 轮询、夹具生成与验证、批量状态文件（`_plan.md`/`_memory.md`）维护、跨文件 grep 排查。

> 唯一例外 (主 Agent 可直接做)：写 ADR/brief、短文件的元数据/状态编辑、输出有界的一次性 `ls`/`git log` 确认——**前提是读入成本极低，超出即委派。**

## 二、三 agent 分工（代码规范在各 agent 定义内，此处只给路由）

> **不在此重复代码规范**——coder 的红→绿测试纪律、四态验收、size 门槛、内存纪律执行，code-reviewer 的红线/审阅维度/回报契约，operator 的四态验收/清单门槛，全部以各自 agent 定义为准：
> `.claude/agents/coder.md`、`.claude/agents/code-reviewer.md`、`.claude/agents/operator.md`。

| 工作 | 委派给 | 本项目典型场景 |
|---|---|---|
| 写/改 **代码**（`src/scrna_integration/*.py`）、写/改 **notebook**（`notebooks/*.ipynb`）、写测试、调试、重构 | **coder** | read_with_manifest / scorers / 各 stage notebook / 夹具脚本 `scripts/make_test_subset.py` |
| **PR diff 审查**（独立会话，与 coder 不同上下文）、verdict + issue 清单 + 回报契约 | **code-reviewer** | 每个 PR 合并前必经；产物是 src/notebook/CI/schema 的 diff 才归它 |
| **非代码执行**：`gh pr merge` 收尾、CI 轮询、worktree 清理、远程分支删除、夹具生成/验证、`_plan.md`/`_memory.md` 批量状态维护、跨文件搜索整理、CLI 调用 | **operator** | PR 合并收尾（见 repo-loop）、跑 `make_test_subset.py` 生成夹具、状态记账 |

**判定**：产物是代码/notebook → coder；产物是 diff 审查 verdict → code-reviewer；产物是被处理后的文件/git 状态/数据 → operator。

## 三、每个 PR 的标准循环（主 Agent 编排，不亲自执行）

1. 主 Agent 按 `_plan.md` 写 **coder brief** → 委派 **coder**（worktree 隔离）实现 + 自验证
2. 主 Agent 读 coder 回报 → 委派 **code-reviewer**（**独立会话**，绝不复用 coder 上下文）审 diff
3. 主 Agent 读 reviewer verdict → 路由决策：
   - `approve` + CI 绿 + 非升级点 → 主 Agent 决定合并，**委派 operator 执行**「等 CI → squash 合并 → 清 worktree → 删本地+**远程**分支 → 更新 `_plan`/`_memory`」（repo-loop 模式，见 `.claude/skills/repo-loop/SKILL.md`）
   - `request-changes` → 主 Agent 据回报契约写再修 brief → 回 coder（新会话）
   - `block`/红线/升级点 → 停下报告 PI
4. **operator 收尾后主 Agent 必须实地核验**（远程分支是否真删、main 是否真合并、状态文件是否真更新）——不只信回报。

**升级 PI 的点**（不自决）：命中红线 / 改保护字段 / 新付费依赖 / 与 `_plan` 方向冲突 / reviewer block / 涉及数据隐私 / 第一波生物学结果（PI 要亲眼看）。

**`git commit` / `git push` 闸门**（主 Agent 每次 commit/push 前自检，缺一步即停止）：
1. [ ] coder 改动已完成并自验证（有回报）
2. [ ] **独立** code-reviewer 已审查（新会话，不复用 coder 上下文），verdict 为 `approve`
3. [ ] coder 的自我评分**不是**独立审查——`QUALITY SCORES` 不能替代 reviewer verdict
4. [ ] reviewer 提出的 issue（P0/P1）已全部修复 / 确认假阳性
5. [ ] 以上 4 条全部满足，才能 `git commit` / `git push`

> 违反纪律案例（2026-06-17）：coder 自评分 0.83 + "PASS" → 主 Agent 跳过独立 reviewer 直接 commit/push。该闸门防止此模式复现。

## 四、本项目特有铁律

- **notebook PR 防超时**：coder **禁止**在单次 turn 内跑 nbconvert 端到端（scVI 训练 + LLM 调用必超时，2026-06-08 实测两次）。改为：静态构建 notebook（nbformat 增量）+ 清 output（`--clear-output`），**运行验证交 PI 在 jupyter 手动跑或 CI**。
- **远程分支清理**：用 `gh api -X DELETE repos/{owner}/{repo}/git/refs/heads/{branch}`（走 HTTPS），比 `git push origin --delete`（SSH）抗网络抖动。`gh pr merge --delete-branch` 只删本地，远程需另删。
- **数据零进 git**：原始数据在 `~/Works/GCPL_scRNA/`（只读），夹具在本地 `data/_subset/`（gitignore，~443M）。只有抽样脚本 + manifest 进 git。worktree 用软链复用主树 `data/_subset/`，加进 `.git/info/exclude`。
- **conda 环境隔离**：所有装包在专用 `scrna-integration` / `scrna-integration-r`，绝不动 base。多 worktree 共享环境用 `PYTHONPATH=src pytest`，禁止 `pip install -e .`（editable 会互相覆盖）。
- **注释中文 + 充分讲 why**（ADR-0009，面向 PI + 非计算机专业学生）。
- **薄框架 + 教学透明**（ADR-0001/0003/0004/0009）：新增任何 helper/抽象，判据是"非 CS 学生打开 notebook 能否逐行看懂"，不只是"填了 scanpy 空白"。
- **跨平台一致性**（ADR-0010）：项目同时跑在 Mac（osx-arm64）+ Linux 服务器（linux-64）。① 版本一致靠 **精确 pin 源 spec + env_parity 诊断脚本 + 人工对齐**（源 spec 用 `==` pin 作为期望基准，两机 `conda env create -f environment.yml` 直接安装，`scripts/env_parity.py snapshot` 留快照，`compare` 出差异供人工决定对齐）② 无法对齐的极少数包登记 `docs/cross-platform-exceptions.md`（异常非常态）③ OS 检测**只能**写在 `src/scrna_integration/platform.py`，notebook/src 其他位置禁出现 `sys.platform`/`os.uname`/平台绝对路径（`/Users/`、`/home/`）④ **conda 环境目录永不进 Syncthing/git**，只同步源 spec+platform.py+异常表+快照 JSON。reviewer 卡这四条。

## 五、关键文件指针

- `_project.md` 元数据 / `_plan.md` PR 计划与状态 / `_memory.md` 项目记忆（含会话保存点，续跑先读它）
- `SPEC.md` 实现规格 / `CONTEXT.md` 术语表 / `docs/adr/` 架构决策（0001-0009）
- `.env.example` LLM provider key+url 模板（06 mLLMCelltype 多模型直连）
