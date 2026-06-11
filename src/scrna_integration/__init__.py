"""scrna-integration-framework: 跨疾病多来源单细胞整合分析框架。

框架代码面积极小——三个公开函数 + 可直接调用的指标函数模块。
一切以面向非计算机专业 PI/学生、能逐行看懂为准（ADR-0009）。
"""

__version__ = "0.0.1"

from scrna_integration.io import inject_genomic_positions as inject_genomic_positions
from scrna_integration.io import read_with_manifest as read_with_manifest
from scrna_integration.markers import load_markers as load_markers
from scrna_integration.platform import check_r_available as check_r_available
from scrna_integration.platform import rscript_bin as rscript_bin
# scorers 模块可直接导入使用：
#   from scrna_integration.scorers import integration_metrics
#   from scrna_integration.scorers import clustering_metrics
# 注意：框架不再 re-export scorers 函数——请在 notebook 中按需从 scorers 模块直接导入。
# platform 模块提供跨平台路径解析（ADR-0010）：
#   from scrna_integration import rscript_bin
#   from scrna_integration.platform import rscript_bin  # 两种写法均可
