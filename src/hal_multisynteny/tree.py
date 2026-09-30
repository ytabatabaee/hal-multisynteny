"""Rooted binary Newick parsing and HAL node-map validation."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path


class TreeError(ValueError):
    """Raised for malformed or unsupported guide trees."""


@dataclass(frozen=True, slots=True)
class TreeNode:
    node_id: str
    label: str | None
    children: tuple[TreeNode, ...] = ()

    @property
    def is_leaf(self) -> bool:
        return not self.children


@dataclass(frozen=True, slots=True)
class TraversalStep:
    node_id: str
    hal_genome: str
    left_child: str
    right_child: str
    left_hal_genome: str
    right_hal_genome: str


@dataclass(frozen=True, slots=True)
class TraversalPlan:
    root: TreeNode
    normalized_newick: str
    tree_sha256: str
    node_map_sha256: str
    node_to_hal: dict[str, str]
    steps: tuple[TraversalStep, ...]


class _Parser:
    def __init__(self, text: str) -> None:
        self.text = text.strip()
        self.index = 0
        self.internal_index = 0

    def parse(self) -> TreeNode:
        if not self.text.endswith(";"):
            raise TreeError("malformed Newick: tree must end with ';'")
        node = self._subtree()
        self._skip_ws()
        if self.index != len(self.text) - 1:
            raise TreeError("malformed Newick: trailing content after root")
        labels: set[str] = set()
        for leaf in _leaves(node):
            if leaf.label in labels:
                raise TreeError(f"duplicate leaf name in guide tree: {leaf.label}")
            labels.add(leaf.label or "")
        return node

    def _subtree(self) -> TreeNode:
        self._skip_ws()
        if self._peek() == "(":
            self.index += 1
            children = [self._subtree()]
            while True:
                self._skip_ws()
                char = self._peek()
                if char == ",":
                    self.index += 1
                    children.append(self._subtree())
                    continue
                if char == ")":
                    self.index += 1
                    break
                raise TreeError("malformed Newick: expected ',' or ')'")
            if len(children) != 2:
                raise TreeError("polytomies are not supported without an explicit resolution policy")
            label = self._label_or_none()
            node_id = label or self._next_internal_id()
            return TreeNode(node_id, label, tuple(children))
        label = self._label_or_none()
        if not label:
            raise TreeError("malformed Newick: leaf nodes must be named")
        return TreeNode(label, label)

    def _label_or_none(self) -> str | None:
        self._skip_ws()
        start = self.index
        while self.index < len(self.text) and self.text[self.index] not in ",():;":
            self.index += 1
        label = self.text[start : self.index].strip()
        if self._peek() == ":":
            self.index += 1
            while self.index < len(self.text) and self.text[self.index] not in ",();":
                self.index += 1
        return label or None

    def _next_internal_id(self) -> str:
        self.internal_index += 1
        return f"internal_{self.internal_index:04d}"

    def _peek(self) -> str:
        return self.text[self.index] if self.index < len(self.text) else ""

    def _skip_ws(self) -> None:
        while self.index < len(self.text) and self.text[self.index].isspace():
            self.index += 1


def parse_newick(text: str) -> TreeNode:
    return _Parser(text).parse()


def _leaves(node: TreeNode) -> tuple[TreeNode, ...]:
    if node.is_leaf:
        return (node,)
    return tuple(leaf for child in node.children for leaf in _leaves(child))


def _postorder(node: TreeNode) -> tuple[TreeNode, ...]:
    return tuple(child_node for child in node.children for child_node in _postorder(child)) + (
        node,
    )


def normalize_newick(node: TreeNode) -> str:
    if node.is_leaf:
        return node.node_id
    return "(" + ",".join(normalize_newick(child) for child in node.children) + ")" + node.node_id


def read_node_map(path: str | Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for line_number, raw in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), start=1):
        if not raw.strip() or raw.startswith("#"):
            continue
        parts = raw.split("\t")
        if len(parts) != 2:
            raise TreeError(f"invalid node map at line {line_number}: expected node_id<TAB>hal_genome")
        node_id, hal_genome = parts
        if not node_id or not hal_genome:
            raise TreeError(f"invalid node map at line {line_number}: empty field")
        if node_id in mapping:
            raise TreeError(f"duplicate node-map entry for {node_id!r}")
        mapping[node_id] = hal_genome
    return mapping


def build_traversal_plan(
    tree_path: str | Path,
    node_map_path: str | Path,
    hal_parent_map: dict[str, str] | None = None,
) -> TraversalPlan:
    tree_text = Path(tree_path).read_text(encoding="utf-8")
    map_text = Path(node_map_path).read_text(encoding="utf-8")
    root = parse_newick(tree_text)
    node_to_hal = read_node_map(node_map_path)
    all_nodes = {node.node_id for node in _postorder(root)}
    missing = sorted(all_nodes - set(node_to_hal))
    extra = sorted(set(node_to_hal) - all_nodes)
    if missing:
        raise TreeError(f"node map is missing guide-tree nodes: {', '.join(missing)}")
    if extra:
        raise TreeError(f"node map contains unknown guide-tree nodes: {', '.join(extra)}")

    steps: list[TraversalStep] = []
    for node in _postorder(root):
        if node.is_leaf:
            continue
        left, right = node.children
        node_hal = node_to_hal[node.node_id]
        left_hal = node_to_hal[left.node_id]
        right_hal = node_to_hal[right.node_id]
        if hal_parent_map is not None:
            for child_id, child_hal in ((left.node_id, left_hal), (right.node_id, right_hal)):
                found = hal_parent_map.get(child_hal)
                if found != node_hal:
                    raise TreeError(
                        "guide-tree/HAL ancestry mismatch: "
                        f"{child_id} maps to HAL genome {child_hal}, parent is {found!r}, "
                        f"expected {node_hal!r}"
                    )
        steps.append(TraversalStep(node.node_id, node_hal, left.node_id, right.node_id, left_hal, right_hal))

    normalized = normalize_newick(root) + ";"
    return TraversalPlan(
        root,
        normalized,
        sha256(normalized.encode("utf-8")).hexdigest(),
        sha256(map_text.encode("utf-8")).hexdigest(),
        node_to_hal,
        tuple(steps),
    )
