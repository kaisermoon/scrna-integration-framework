"""scrna-integration-framework: cross-disease, multi-source scRNA-seq integration."""

__version__ = "0.0.1"

from scrna_integration.io import read_with_manifest as read_with_manifest
from scrna_integration.markers import load_markers as load_markers
from scrna_integration.sweep import sweep as sweep
