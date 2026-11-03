import sys
import argparse
from collections import defaultdict


class Interval:
    def __init__(self, low=None, high=None):
        self.low = low
        self.high = high

    def empty(self):
        return self.low is None and self.high is None

    def copy(self):
        return Interval(self.low, self.high)

    def conflicting(self, b, state):
        # returns true if this interval conflicts with back-edge b
        return not self.empty() and state.lowpt[self.high] > state.lowpt[b]

    def __repr__(self):
        return f"[{self.low}, {self.high}]" if not self.empty() else "[]"


class ConflictPair:
    def __init__(self, left=None, right=None):
        self.left = left if left is not None else Interval()
        self.right = right if right is not None else Interval()

    def swap(self):
        self.left, self.right = self.right, self.left

    def lowest(self, state):
        if self.left.empty():
            return state.lowpt[self.right.low]
        if self.right.empty():
            return state.lowpt[self.left.low]
        return min(state.lowpt[self.left.low], state.lowpt[self.right.low])

    def __repr__(self):
        return f"ConflictPair(L={self.left}, R={self.right})"


def top_of_stack(stack):
    return stack[-1] if stack else None


class LRPlanarityTest:
    def __init__(self, vertices, edges, verbose=False):
        self.raw_vertices = list(vertices)
        self.raw_edges = list(edges)
        self.verbose = verbose

        # map vertex names/labels to numbers 0..V-1
        self.node_to_id = {node: i for i, node in enumerate(self.raw_vertices)}
        self.id_to_node = {i: node for node, i in self.node_to_id.items()}
        self.V = len(self.raw_vertices)

        self.adj = defaultdict(list)
        self.clean_edges = []
        seen = set()

        self.self_loops = 0
        self.parallel_edges = 0

        # remove self-loops and duplicate edges
        for u, v in self.raw_edges:
            if u not in self.node_to_id:
                self.node_to_id[u] = len(self.node_to_id)
                self.id_to_node[self.node_to_id[u]] = u
                self.V += 1
            if v not in self.node_to_id:
                self.node_to_id[v] = len(self.node_to_id)
                self.id_to_node[self.node_to_id[v]] = v
                self.V += 1

            uid = self.node_to_id[u]
            vid = self.node_to_id[v]

            if uid == vid:
                self.self_loops += 1
                continue

            edge = (min(uid, vid), max(uid, vid))
            if edge in seen:
                self.parallel_edges += 1
                continue

            seen.add(edge)
            self.clean_edges.append(edge)
            self.adj[uid].append(vid)
            self.adj[vid].append(uid)

        self.E = len(self.clean_edges)

        # dfs tree and lowpoints
        self.roots = []
        self.height = defaultdict(lambda: None)
        self.lowpt = {}
        self.lowpt2 = {}
        self.nesting_depth = {}
        self.parent_edge = defaultdict(lambda: None)

        # oriented edges
        self.DG_edges = set()
        self.DG_adj = defaultdict(list)
        self.ordered_adjs = {}

        # conflict resolution structures
        self.ref = defaultdict(lambda: None)
        self.side = defaultdict(lambda: 1)  # 1 = right, -1 = left
        self.S = []
        self.stack_bottom = {}
        self.lowpt_edge = {}

    def run(self):
        if self.verbose:
            print(f"Testing graph: V={self.V}, E={self.E}")
            if self.self_loops or self.parallel_edges:
                print(f"  (Ignored {self.self_loops} self-loops and {self.parallel_edges} parallel edges)")

        # Euler's formula check: planar simple graph has E <= 3V - 6
        if self.V >= 3 and self.E > 3 * self.V - 6:
            reason = f"Euler formula check failed: E ({self.E}) > 3V - 6 ({3 * self.V - 6})"
            return False, reason, {"V": self.V, "E": self.E}

        # orient graph using DFS
        for v in range(self.V):
            if self.height[v] is None:
                self.height[v] = 0
                self.roots.append(v)
                self.dfs_orientation(v)

        # sort adjacency list by nesting depth
        for v in range(self.V):
            self.ordered_adjs[v] = sorted(
                self.DG_adj[v], key=lambda x: self.nesting_depth[(v, x)]
            )

        # test planarity on each DFS component
        for r in self.roots:
            if not self.dfs_testing(r):
                reason = "Conflict detected: back-edges cannot be partitioned without crossing"
                return False, reason, {"V": self.V, "E": self.E}

        # resolve final Left/Right orientations
        for e in list(self.DG_edges):
            self.sign_recursive(e)

        lr_partition = {}
        for (u, v) in self.DG_edges:
            orig_u = self.id_to_node[u]
            orig_v = self.id_to_node[v]
            side_str = "Right" if self.side[(u, v)] == 1 else "Left"
            is_tree = (self.parent_edge[v] == (u, v))
            lr_partition[(orig_u, orig_v)] = {
                "type": "Tree Edge" if is_tree else "Back Edge",
                "side": side_str
            }

        return True, "Graph is planar", {"lr_partition": lr_partition}

    def dfs_orientation(self, v):
        # iterative dfs to avoid recursion limit
        dfs_stack = [v]
        ind = defaultdict(lambda: 0)
        skip_init = defaultdict(lambda: False)

        while dfs_stack:
            curr = dfs_stack.pop()
            e = self.parent_edge[curr]

            for w in self.adj[curr][ind[curr]:]:
                vw = (curr, w)
                wv = (w, curr)

                if not skip_init[vw]:
                    if vw in self.DG_edges or wv in self.DG_edges:
                        ind[curr] += 1
                        continue

                    self.DG_edges.add(vw)
                    self.DG_adj[curr].append(w)

                    self.lowpt[vw] = self.height[curr]
                    self.lowpt2[vw] = self.height[curr]

                    if self.height[w] is None:
                        # tree edge
                        self.parent_edge[w] = vw
                        self.height[w] = self.height[curr] + 1
                        dfs_stack.append(curr)
                        dfs_stack.append(w)
                        skip_init[vw] = True
                        break
                    else:
                        # back edge
                        self.lowpt[vw] = self.height[w]

                self.nesting_depth[vw] = 2 * self.lowpt[vw]
                if self.lowpt2[vw] < self.height[curr]:
                    self.nesting_depth[vw] += 1

                if e is not None:
                    if self.lowpt[vw] < self.lowpt[e]:
                        self.lowpt2[e] = min(self.lowpt[e], self.lowpt2[vw])
                        self.lowpt[e] = self.lowpt[vw]
                    elif self.lowpt[vw] > self.lowpt[e]:
                        self.lowpt2[e] = min(self.lowpt2[e], self.lowpt[vw])
                    else:
                        self.lowpt2[e] = min(self.lowpt2[e], self.lowpt2[vw])

                ind[curr] += 1

    def dfs_testing(self, v):
        dfs_stack = [v]
        ind = defaultdict(lambda: 0)
        skip_init = defaultdict(lambda: False)

        while dfs_stack:
            curr = dfs_stack.pop()
            e = self.parent_edge[curr]
            skip_final = False

            for w in self.ordered_adjs[curr][ind[curr]:]:
                ei = (curr, w)

                if not skip_init[ei]:
                    self.stack_bottom[ei] = top_of_stack(self.S)

                    if ei == self.parent_edge[w]:
                        # tree edge
                        dfs_stack.append(curr)
                        dfs_stack.append(w)
                        skip_init[ei] = True
                        skip_final = True
                        break
                    else:
                        # back edge
                        self.lowpt_edge[ei] = ei
                        self.S.append(ConflictPair(right=Interval(ei, ei)))

                if self.lowpt[ei] < self.height[curr]:
                    if w == self.ordered_adjs[curr][0]:
                        self.lowpt_edge[e] = self.lowpt_edge[ei]
                    else:
                        if not self.add_constraints(ei, e):
                            return False

                ind[curr] += 1

            if not skip_final:
                if e is not None:
                    self.remove_back_edges(e)

        return True

    def add_constraints(self, ei, e):
        P = ConflictPair()

        # merge return edges of ei into P.right
        while True:
            Q = self.S.pop()
            if not Q.left.empty():
                Q.swap()
            if not Q.left.empty():
                return False

            if self.lowpt[Q.right.low] > self.lowpt[e]:
                if P.right.empty():
                    P.right = Q.right.copy()
                else:
                    self.ref[P.right.low] = Q.right.high
                P.right.low = Q.right.low
            else:
                self.ref[Q.right.low] = self.lowpt_edge[e]

            if top_of_stack(self.S) == self.stack_bottom[ei]:
                break

        # merge conflicting return edges into P.left
        while (top_of_stack(self.S) and (
            top_of_stack(self.S).left.conflicting(ei, self) or
            top_of_stack(self.S).right.conflicting(ei, self)
        )):
            Q = self.S.pop()
            if Q.right.conflicting(ei, self):
                Q.swap()
            if Q.right.conflicting(ei, self):
                return False

            self.ref[P.right.low] = Q.right.high
            if Q.right.low is not None:
                P.right.low = Q.right.low

            if P.left.empty():
                P.left = Q.left.copy()
            else:
                self.ref[P.left.low] = Q.left.high
            P.left.low = Q.left.low

        if not (P.left.empty() and P.right.empty()):
            self.S.append(P)

        return True

    def remove_back_edges(self, e):
        u = e[0]

        # trim back edges returning to u
        while self.S and top_of_stack(self.S).lowest(self) == self.height[u]:
            P = self.S.pop()
            if P.left.low is not None:
                self.side[P.left.low] = -1

        if self.S:
            P = self.S.pop()
            while P.left.high is not None and P.left.high[1] == u:
                P.left.high = self.ref[P.left.high]
            if P.left.high is None and P.left.low is not None:
                self.ref[P.left.low] = P.right.low
                self.side[P.left.low] = -1
                P.left.low = None

            while P.right.high is not None and P.right.high[1] == u:
                P.right.high = self.ref[P.right.high]
            if P.right.high is None and P.right.low is not None:
                self.ref[P.right.low] = P.left.low
                self.side[P.right.low] = -1
                P.right.low = None

            self.S.append(P)

        if self.lowpt[e] < self.height[u]:
            hl = top_of_stack(self.S).left.high
            hr = top_of_stack(self.S).right.high

            if hl is not None and (hr is None or self.lowpt[hl] > self.lowpt[hr]):
                self.ref[e] = hl
            else:
                self.ref[e] = hr

    def sign_recursive(self, e):
        # resolve relative orientation to absolute side
        if self.ref[e] is not None:
            self.side[e] = self.side[e] * self.sign_recursive(self.ref[e])
            self.ref[e] = None
        return self.side[e]


def parse_graphs_from_file(filepath):
    graphs = []
    with open(filepath, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]

    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("GRAPH"):
            parts = line.split(maxsplit=1)
            name = parts[1] if len(parts) > 1 else f"Graph_{len(graphs)+1}"
            i += 1
            meta = lines[i].split()
            v_count = int(meta[0])
            e_count = int(meta[1])
            i += 1
            edges = []
            vertices = set()
            for _ in range(e_count):
                if i < len(lines) and lines[i] != "END":
                    u, v = lines[i].split()[:2]
                    edges.append((u, v))
                    vertices.add(u)
                    vertices.add(v)
                    i += 1
            if i < len(lines) and lines[i] == "END":
                i += 1
            all_verts = list(vertices)
            if len(all_verts) < v_count:
                for idx in range(v_count):
                    all_verts.append(str(idx))
                all_verts = list(dict.fromkeys(all_verts))[:v_count]
            graphs.append((name, all_verts, edges))
        else:
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
                v_count = int(parts[0])
                e_count = int(parts[1])
                name = f"Graph_{len(graphs)+1}"
                i += 1
                edges = []
                vertices = set()
                for _ in range(e_count):
                    if i < len(lines):
                        ep = lines[i].split()
                        if len(ep) >= 2:
                            u, v = ep[0], ep[1]
                            edges.append((u, v))
                            vertices.add(u)
                            vertices.add(v)
                        i += 1
                all_verts = list(vertices)
                if len(all_verts) < v_count:
                    for idx in range(v_count):
                        all_verts.append(str(idx))
                    all_verts = list(dict.fromkeys(all_verts))[:v_count]
                graphs.append((name, all_verts, edges))
            else:
                i += 1

    return graphs


def run_demo():
    demo_cases = [
        ("Complete Graph K4", 4, [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (1, 3)]),
        ("Complete Graph K5", 5, [
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
            (0, 2), (2, 4), (4, 1), (1, 3), (3, 0)
        ]),
        ("Utility Graph K3,3", 6, [
            (0, 3), (0, 4), (0, 5),
            (1, 3), (1, 4), (1, 5),
            (2, 3), (2, 4), (2, 5)
        ]),
        ("3D Cube Graph (Q3)", 8, [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7)
        ]),
        ("Petersen Graph", 10, [
            (0, 1), (0, 4), (0, 5), (1, 2), (1, 6),
            (2, 3), (2, 7), (3, 4), (3, 8), (4, 9),
            (5, 7), (5, 8), (6, 8), (6, 9), (7, 9)
        ]),
        ("Pentagonal Prism (C5 x K2)", 10, [
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
            (5, 6), (6, 7), (7, 8), (8, 9), (9, 5),
            (0, 5), (1, 6), (2, 7), (3, 8), (4, 9)
        ]),
        ("W9 Lecture Graph (Slide 4)", 7, [
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 0),
            (1, 6), (2, 6), (2, 5), (3, 5)
        ]),
    ]

    for name, v_count, edges in demo_cases:
        print(f"\n>> Evaluating: {name} (|V|={v_count}, |E|={len(edges)})")
        tester = LRPlanarityTest(list(range(v_count)), edges, verbose=False)
        is_planar, reason, details = tester.run()

        status_str = "[ PLANAR ]" if is_planar else "[ NON-PLANAR ]"
        print(f"   Status : {status_str}")
        print(f"   Reason : {reason}")
        if is_planar and "lr_partition" in details:
            print("   LR-Partition of Edges:")
            for e, info in list(details["lr_partition"].items())[:6]:
                print(f"     - Edge {e[0]} -> {e[1]}: {info['type']:10} (Side: {info['side']})")
            if len(details["lr_partition"]) > 6:
                print(f"     - ... ({len(details['lr_partition']) - 6} more edges)")


def main():
    parser = argparse.ArgumentParser(description="Left-Right Planarity Test")
    parser.add_argument("input_file", nargs="?", help="Path to graph input file")
    parser.add_argument("--demo", action="store_true", help="Run demo graphs")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    args = parser.parse_args()

    if args.demo or not args.input_file:
        run_demo()
        if not args.input_file:
            print("\nTip: To test with a file, run:")
            print("  python planarity_lr.py input_samples.txt [-v]")
            return

    if args.input_file:
        graphs = parse_graphs_from_file(args.input_file)
        if not graphs:
            print("Error: No valid graphs found in file.")
            sys.exit(1)

        print(f"Found {len(graphs)} graph(s) to evaluate.\n" + "-" * 50)
        for name, vertices, edges in graphs:
            print(f"\n>> Graph: '{name}' (V={len(vertices)}, E={len(edges)})")
            tester = LRPlanarityTest(vertices, edges, verbose=args.verbose)
            is_planar, reason, details = tester.run()

            tag = "[ PLANAR ]" if is_planar else "[ NON-PLANAR ]"
            print(f"   Verdict: {tag}")
            print(f"   Details: {reason}")
            if is_planar and "lr_partition" in details:
                print(f"   LR-Partition ({len(details['lr_partition'])} edges):")
                for e, info in details["lr_partition"].items():
                    print(f"     * {e[0]} -> {e[1]}: {info['type']} ({info['side']})")


if __name__ == "__main__":
    main()
