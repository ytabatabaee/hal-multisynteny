"""Ancestry-aware multi-genome synteny block construction."""

from .builder import BuildConfig, BuildResult, build_blocks
from .models import AlignmentRun, AncestralBlock, BlockOccurrence, ParentMappedRun
from .reconcile import ReconcileConfig, ReconcileResult, reconcile_node

__all__ = [
    "AlignmentRun",
    "AncestralBlock",
    "BlockOccurrence",
    "BuildConfig",
    "BuildResult",
    "ParentMappedRun",
    "ReconcileConfig",
    "ReconcileResult",
    "build_blocks",
    "reconcile_node",
]

__version__ = "0.2.0"
