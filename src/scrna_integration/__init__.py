"""scrna-integration-framework: 跨疾病多来源单细胞整合分析框架。

src/notebook 边界铁律（2026-07-10 PI 定稿）：进 src/ 的仅为技术管道。
框架代码面积极小——面向非计算机专业 PI/学生、能逐行看懂为准（ADR-0009）。
"""

__version__ = "0.0.1"

from scrna_integration.io import inject_genomic_positions as inject_genomic_positions
from scrna_integration.io import sync_gene_ids as sync_gene_ids

# LLM 统一调用（批 5，llm_config.py）
from scrna_integration.llm_config import apply_mllmcelltype_patches as apply_mllmcelltype_patches
from scrna_integration.llm_config import build_mllmcelltype_config as build_mllmcelltype_config
from scrna_integration.llm_config import call_llm_for_annotation as call_llm_for_annotation
from scrna_integration.llm_config import (
    extract_json_from_llm_response as extract_json_from_llm_response,
)

# platform 模块提供跨平台路径解析（ADR-0010）：
#   from scrna_integration import rscript_bin
#   from scrna_integration.platform import rscript_bin  # 两种写法均可
from scrna_integration.platform import check_r_available as check_r_available
from scrna_integration.platform import rscript_bin as rscript_bin

# 运行契约与状态机（P0 契约地基）
from scrna_integration.run_contract import MethodStatus as MethodStatus
from scrna_integration.run_contract import StageStatus as StageStatus
from scrna_integration.run_contract import aggregate_method_status as aggregate_method_status
from scrna_integration.run_contract import determine_stage_status as determine_stage_status
from scrna_integration.run_contract import (
    validate_expression_contract as validate_expression_contract,
)
