"""Ancestry-aware multi-genome synteny block construction."""

from .builder import BuildConfig, BuildResult, build_blocks
from .models import AlignmentRun, AncestralBlock, BlockOccurrence

__all__ = [
    "AlignmentRun",
    "AncestralBlock",
    "BlockOccurrence",
    "BuildConfig",
    "BuildResult",
    "build_blocks",
]

__version__ = "0.1.0"
