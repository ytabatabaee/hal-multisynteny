"""Ancestry-aware multi-genome synteny block construction."""

from .builder import BuildConfig, BuildResult, build_blocks
from .models import (
    AlignmentRun,
    AncestralBlock,
    BlockOccurrence,
    EdgeMappingRun,
    LeafOccurrence,
    NodeBlock,
    NodeBlockOccurrence,
    ParentMappedRun,
    UnmappedEdgeEvidence,
)
from .reconcile import ReconcileConfig, ReconcileResult, reconcile_node

__all__ = [
    "AlignmentRun",
    "AncestralBlock",
    "BlockOccurrence",
    "BuildConfig",
    "BuildResult",
    "EdgeMappingRun",
    "LeafOccurrence",
    "NodeBlock",
    "NodeBlockOccurrence",
    "ParentMappedRun",
    "ReconcileConfig",
    "ReconcileResult",
    "UnmappedEdgeEvidence",
    "build_blocks",
    "reconcile_node",
]

__version__ = "0.3.1"
