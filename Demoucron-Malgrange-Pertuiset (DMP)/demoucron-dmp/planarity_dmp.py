import sys
import argparse
from collections import defaultdict, deque


class DMPPlanarityTest:
    """
    Planarity testing using the Demoucron-Malgrange-Pertuiset (DMP) algorithm (1964).
    
    The algorithm is an incremental cycle/path addition method:
    1. Preprocessing: Remove self-loops and parallel edges.
    2. Quick Euler formula bound check: E <= 3V - 6 for V >= 3.
    3. Decompose graph into 2-connected components (blocks).
    4. For each non-trivial block:
       a. Find an initial simple cycle C. C defines two initial faces (interior and exterior).
       b. Repeatedly identify 'bridges' (fragments) of un-embedded edges.
       c. Compute admissible faces for each bridge (faces containing all bridge attachments).
       d. If any bridge has 0 admissible faces, the graph is NON-PLANAR.
       e. Select a bridge with the fewest admissible faces (forcing constraint).
       f. Find an alpha-path within the bridge connecting two attachments.
       g. Embed the alpha-path into one of its admissible faces, splitting it into two new faces.
    5. If all blocks are successfully embedded, the graph is PLANAR.
    """

    def __init__(self, vertices, edges, verbose=False):
        self.raw_vertices = list(vertices)
        self.raw_edges = list(edges)
        self.verbose = verbose

        # Map vertex identifiers/labels to indices 0..V-1
        self.node_to_id = {node: i for i, node in enumerate(self.raw_vertices)}
        self.id_to_node = {i: node for node, i in self.node_to_id.items()}
        self.V = len(self.raw_vertices)

        self.clean_edges = []
        seen = set()
        self.self_loops = 0
        self.parallel_edges = 0

        # Clean multi-edges and self-loops
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

    def _log(self, msg):
        if self.verbose:
            print(f"  [DMP] {msg}")

    def get_biconnected_components(self):
        """
        Decomposes the simplified graph into biconnected components (blocks)
        using Tarjan's bridge/articulation DFS algorithm.
        """
        adj = defaultdict(list)
        for u, v in self.clean_edges:
            adj[u].append(v)
            adj[v].append(u)

        visited = set()
        tin = {}
        low = {}
        timer = 0
        stack = []
        blocks = []

        def dfs(u, p=-1):
            nonlocal timer
            visited.add(u)
            tin[u] = low[u] = timer
            timer += 1

            for v in adj[u]:
                if v == p:
                    continue
                if v in visited:
                    if tin[v] < tin[u]:
                        stack.append((u, v))
                        low[u] = min(low[u], tin[v])
                else:
                    stack.append((u, v))
                    dfs(v, u)
                    low[u] = min(low[u], low[v])
                    if low[v] >= tin[u]:
                        # u is an articulation point (or DFS tree root)
                        block = []
                        while stack:
                            edge = stack.pop()
                            block.append(edge)
                            if edge == (u, v):
                                break
                        blocks.append(block)

        for i in range(self.V):
            if i not in visited:
                dfs(i)

        return blocks

    def _find_initial_cycle(self, block_nodes, adj):
        """
        Finds a simple cycle of length >= 3 in a 2-connected block.
        """
        start = next(iter(block_nodes))
        parent = {start: None}
        cycle = []

        def dfs(u, p):
            nonlocal cycle
            if cycle:
                return
            for v in adj[u]:
                if v == p:
                    continue
                if v in parent:
                    # Back-edge found -> cycle detected!
                    path = [u]
                    curr = u
                    while curr != v:
                        curr = parent[curr]
                        path.append(curr)
                    cycle = path
                    return
                else:
                    parent[v] = u
                    dfs(v, u)
                    if cycle:
                        return

        dfs(start, None)
        return cycle

    def _test_block_dmp(self, block_edges):
        """
        Executes Demoucron-Malgrange-Pertuiset on a single 2-connected block.
        Returns (is_planar, reason, faces).
        """
        if len(block_edges) <= 3:
            # Blocks with <= 3 edges are either single edges, self-contained cycles C3, etc.
            return True, "Trivially planar small block", []

        block_nodes = set()
        adj = defaultdict(set)
        for u, v in block_edges:
            block_nodes.add(u)
            block_nodes.add(v)
            adj[u].add(v)
            adj[v].add(u)

        V = len(block_nodes)
        E = len(block_edges)

        # Euler's formula inequality for planar simple graphs: E <= 3V - 6
        if E > 3 * V - 6:
            return False, f"Euler formula check failed for block: E ({E}) > 3V - 6 ({3 * V - 6})", []

        # Find initial simple cycle C
        cycle = self._find_initial_cycle(block_nodes, adj)
        if not cycle or len(cycle) < 3:
            return True, "No cycle of length >= 3 in block", []

        self._log(f"Initial cycle: {' -> '.join(str(self.id_to_node[x]) for x in cycle)} -> {self.id_to_node[cycle[0]]}")

        # H represents currently embedded subgraph
        H_nodes = set(cycle)
        H_edges = set()
        k = len(cycle)
        for i in range(k):
            u, v = cycle[i], cycle[(i + 1) % k]
            H_edges.add((min(u, v), max(u, v)))

        # Faces in planar embedding. Cycle C partitions the sphere/plane into 2 faces:
        # Face 1: clockwise traversal of C
        # Face 2: counter-clockwise traversal of C
        faces = [
            list(cycle),
            list(reversed(cycle))
        ]

        all_block_edges = set((min(u, v), max(u, v)) for u, v in block_edges)

        step = 0
        while len(H_edges) < len(all_block_edges):
            step += 1
            rem_edges = all_block_edges - H_edges

            # Build adjacency of remaining edges
            rem_adj = defaultdict(set)
            for u, v in rem_edges:
                rem_adj[u].add(v)
                rem_adj[v].add(u)

            # Decompose remaining edges into bridges relative to H:
            bridges = []

            # 1. Chords: single edges where both endpoints belong to H
            for u, v in rem_edges:
                if u in H_nodes and v in H_nodes:
                    bridges.append({
                        'type': 'chord',
                        'edges': [(u, v)],
                        'attachments': {u, v},
                        'u': u,
                        'v': v
                    })

            # 2. Connected components in block_nodes \ H_nodes
            unvisited_rem_nodes = (block_nodes - H_nodes) & set(rem_adj.keys())
            visited_comp = set()

            for node in unvisited_rem_nodes:
                if node in visited_comp:
                    continue
                comp = []
                q = deque([node])
                visited_comp.add(node)
                while q:
                    curr = q.popleft()
                    comp.append(curr)
                    for nxt in rem_adj[curr]:
                        if nxt not in H_nodes and nxt not in visited_comp:
                            visited_comp.add(nxt)
                            q.append(nxt)

                # Attachments: vertices in H adjacent to this component
                attachments = set()
                b_edges = []
                for curr in comp:
                    for nxt in rem_adj[curr]:
                        e = (min(curr, nxt), max(curr, nxt))
                        if nxt in H_nodes:
                            attachments.add(nxt)
                            if e in rem_edges:
                                b_edges.append(e)
                        elif curr < nxt:
                            b_edges.append(e)

                bridges.append({
                    'type': 'component',
                    'comp': comp,
                    'edges': b_edges,
                    'attachments': attachments
                })

            # Determine admissible faces for each bridge
            candidate_bridge = None
            min_admissible = float('inf')
            admissible_faces_for_candidate = []

            for br in bridges:
                att = br['attachments']
                adm_faces = []
                for f_idx, f in enumerate(faces):
                    f_nodes = set(f)
                    if att.issubset(f_nodes):
                        adm_faces.append(f_idx)

                if len(adm_faces) == 0:
                    att_labels = [str(self.id_to_node[x]) for x in att]
                    return False, f"Bridge with attachments {att_labels} has 0 admissible faces (cannot be drawn without crossing)", []

                if len(adm_faces) < min_admissible:
                    min_admissible = len(adm_faces)
                    candidate_bridge = br
                    admissible_faces_for_candidate = adm_faces
                    if min_admissible == 1:
                        # Forced choice: bridge has only 1 admissible face
                        break

            # Find alpha-path in candidate bridge connecting two attachments
            br = candidate_bridge
            target_face_idx = admissible_faces_for_candidate[0]
            target_face = faces[target_face_idx]

            if br['type'] == 'chord':
                u, v = br['u'], br['v']
                alpha_path = [u, v]
            else:
                comp_set = set(br['comp'])
                att_list = list(br['attachments'])
                start_att = att_list[0]

                parent = {start_att: None}
                q = deque([start_att])
                end_att = None

                while q:
                    curr = q.popleft()
                    if curr != start_att and curr in br['attachments']:
                        end_att = curr
                        break

                    for nxt in rem_adj[curr]:
                        if (nxt in comp_set or (nxt in br['attachments'] and nxt != start_att)) and nxt not in parent:
                            parent[nxt] = curr
                            q.append(nxt)

                path = []
                curr = end_att
                while curr is not None:
                    path.append(curr)
                    curr = parent[curr]
                alpha_path = list(reversed(path))

            u, v = alpha_path[0], alpha_path[-1]
            self._log(
                f"Step {step}: Embedding alpha-path {' -> '.join(str(self.id_to_node[x]) for x in alpha_path)} "
                f"(Bridge attachments: {len(br['attachments'])}, Admissible faces: {min_admissible})"
            )

            # Split target face into two new faces
            idx_u = target_face.index(u)
            idx_v = target_face.index(v)

            N = len(target_face)
            # Segment along cycle boundary from u to v
            seg_u_to_v = []
            curr = idx_u
            while True:
                seg_u_to_v.append(target_face[curr])
                if curr == idx_v:
                    break
                curr = (curr + 1) % N

            # Segment along cycle boundary from v to u
            seg_v_to_u = []
            curr = idx_v
            while True:
                seg_v_to_u.append(target_face[curr])
                if curr == idx_u:
                    break
                curr = (curr + 1) % N

            internal_path = alpha_path[1:-1]

            # Face 1: u -> v along boundary, then reverse of internal alpha-path back to u
            new_face_1 = seg_u_to_v + list(reversed(internal_path))
            # Face 2: v -> u along boundary, then forward of internal alpha-path to v
            new_face_2 = seg_v_to_u + internal_path

            faces.pop(target_face_idx)
            faces.append(new_face_1)
            faces.append(new_face_2)

            # Update embedded graph H
            for node in alpha_path:
                H_nodes.add(node)
            for i in range(len(alpha_path) - 1):
                e = (min(alpha_path[i], alpha_path[i + 1]), max(alpha_path[i], alpha_path[i + 1]))
                H_edges.add(e)

        return True, "Embedded successfully", faces

    def run(self):
        """
        Executes the full planarity test on the input graph.
        Returns:
            (is_planar: bool, reason: str, details: dict)
        """
        if self.V <= 3:
            return True, "Graph with <= 3 vertices is trivially planar", {"faces": []}

        E = len(self.clean_edges)
        if E > 3 * self.V - 6:
            return False, f"Euler formula check failed: E ({E}) > 3V - 6 ({3 * self.V - 6})", {}

        blocks = self.get_biconnected_components()
        self._log(f"Decomposed graph into {len(blocks)} biconnected component(s).")

        all_faces = []
        for idx, block in enumerate(blocks):
            if not block:
                continue
            is_planar, reason, faces = self._test_block_dmp(block)
            if not is_planar:
                return False, f"Block {idx + 1}/{len(blocks)} non-planar: {reason}", {}
            
            # Map face vertex IDs back to original labels
            for f in faces:
                labeled_face = [self.id_to_node[x] for x in f]
                all_faces.append(labeled_face)

        return True, "Graph is planar (all biconnected blocks successfully embedded)", {"faces": all_faces}


def parse_graphs_from_file(filename):
    """
    Parses multiple graph specifications from an input text file.
    Supports format:
      GRAPH <name>
      <V> <E>
      u v
      ...
      END
    """
    graphs = []
    with open(filename, "r") as f:
        lines = [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]

    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("GRAPH"):
            parts = line.split()
            name = parts[1] if len(parts) > 1 else f"Graph_{len(graphs)+1}"
            i += 1
            if i >= len(lines):
                break
            dim_parts = lines[i].split()
            v_count = int(dim_parts[0])
            e_count = int(dim_parts[1])
            i += 1

            edges = []
            vertices = set()
            while i < len(lines) and lines[i] != "END":
                ep = lines[i].split()
                if len(ep) >= 2:
                    u, v = ep[0], ep[1]
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
        tester = DMPPlanarityTest(list(range(v_count)), edges, verbose=False)
        is_planar, reason, details = tester.run()

        status_str = "[ PLANAR ]" if is_planar else "[ NON-PLANAR ]"
        print(f"   Status : {status_str}")
        print(f"   Reason : {reason}")
        if is_planar and "faces" in details and details["faces"]:
            print(f"   Planar Faces ({len(details['faces'])} faces):")
            for idx, face in enumerate(details["faces"][:5]):
                face_str = " -> ".join(map(str, face)) + f" -> {face[0]}"
                print(f"     * Face {idx+1}: {face_str}")
            if len(details["faces"]) > 5:
                print(f"     * ... ({len(details['faces']) - 5} more faces)")


def main():
    parser = argparse.ArgumentParser(description="Demoucron-Malgrange-Pertuiset (DMP) Planarity Test")
    parser.add_argument("input_file", nargs="?", help="Path to graph input file")
    parser.add_argument("--demo", action="store_true", help="Run demo graphs")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose step-by-step trace")

    args = parser.parse_args()

    if args.demo or not args.input_file:
        run_demo()
        if not args.input_file:
            print("\nTip: To test with a file, run:")
            print("  python planarity_dmp.py input_samples.txt [-v]")
            return

    if args.input_file:
        graphs = parse_graphs_from_file(args.input_file)
        if not graphs:
            print("Error: No valid graphs found in file.")
            sys.exit(1)

        print(f"Found {len(graphs)} graph(s) to evaluate.\n" + "-" * 50)
        for name, vertices, edges in graphs:
            print(f"\n>> Graph: '{name}' (V={len(vertices)}, E={len(edges)})")
            tester = DMPPlanarityTest(vertices, edges, verbose=args.verbose)
            is_planar, reason, details = tester.run()

            tag = "[ PLANAR ]" if is_planar else "[ NON-PLANAR ]"
            print(f"   Verdict: {tag}")
            print(f"   Details: {reason}")
            if is_planar and "faces" in details and details["faces"]:
                print(f"   Planar Faces ({len(details['faces'])} faces):")
                for f_idx, face in enumerate(details["faces"]):
                    face_str = " -> ".join(map(str, face)) + f" -> {face[0]}"
                    print(f"     * Face {f_idx+1}: {face_str}")


if __name__ == "__main__":
    main()
