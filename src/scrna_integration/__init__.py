"""scrna-integration-framework: cross-disease, multi-source scRNA-seq integration."""

__version__ = "0.0.1"

from scrna_integration.markers import load_markers as load_markers
from scrna_integration.sweep import sweep as sweep

# read_with_manifest will be available after PR-1 merges.
# Use try/except so the package remains importable before that.
try:
    from scrna_integration.io import read_with_manifest  # noqa: F401
except ImportError:
    pass
