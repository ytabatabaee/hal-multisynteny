# Evidence-based chaining proposal

The audit shows that Method 3's short blocks are not primarily caused by a lack
of shared evidence. On 10 Mb, 306,908 of 365,420 atoms are
`shared_consistent`, but the current common-refinement policy preserves every
atomic boundary. The largest primary rejection class is the inherited/
classification barrier (83,302 boundaries; 76.2% of rejected joins), followed by
non-monotonic child-coordinate order (25,385 primary boundaries; 23.2%). Raw
secondary evidence includes 34,648 chromosome transitions, 34,628 orientation
transitions, and 61,808/48,684 forward/reverse collinearity failures. These
categories overlap; they must not be added as mutually exclusive counts.

## Recommended hybrid

1. On each direct edge, normalize native segments into ordered anchors. Keep all
   source IDs and mapping alternatives, but merge adjacent components only when
   parent chromosome, child chromosome, relative orientation, and monotonic
   order agree. Compute parent and orientation-aware child gaps explicitly.
2. At each ancestral node, build candidate chains from both child tracks. Use a
   dynamic-programming weighted interval chain rather than greedy adjacent atom
   joins. Score aligned supported length positively and parent/child gaps
   negatively; apply hard barriers for chromosome changes, orientation
   transitions, non-monotonic order, unresolved ambiguity/duplication, and
   inherited structural boundaries.
3. Emit chain-to-atom provenance and preserve all rejected alternatives. The
   existing reconciliation/classification layer then evaluates the selected
   compatible chains without changing biological statuses.

Pairwise preprocessing is cheaper and reusable but can create incompatible
chains independently in the two descendants. Joint node-level chaining sees
both descendants and is safer for multispecies support but costs more and can
repeat work. The hybrid keeps pairwise indexing reusable and performs the final
joint choice at each ancestor.

For N ordered anchors, dynamic programming is O(N) for a single compatible
track after sorting, or O(N log N) with interval predecessor queries. Joint
multi-track candidate enumeration is O(N log N + K) for K compatible overlaps;
retaining alternatives can increase K but does not require a leaf-by-leaf
Cartesian expansion. The current implementation should remain the production
baseline until this proposal is validated on controlled fixtures and compared
against the same HAL evidence.

Gap relaxation alone is insufficient: in the 100-kb test, raising the threshold
from 1 kb to 10 kb changes 1,014 to 1,009 blocks, while classification and
structural barriers remain. Large gap thresholds should therefore remain an
explicit exploratory parameter, never a replacement for anchor ordering and
barrier checks.
