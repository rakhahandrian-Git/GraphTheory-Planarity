import sys
from planarity_lr import LRPlanarityTest


def run_tests():
    test_cases = [
        ("Single Vertex", 1, [], True),
        ("Single Edge", 2, [(0, 1)], True),
        ("Isolated Vertices (V=4)", 4, [], True),
        ("Triangle Cycle C3", 3, [(0, 1), (1, 2), (2, 0)], True),
        ("Cycle C5", 5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)], True),
        ("Complete Graph K4", 4, [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (1, 3)], True),
        ("Complete Graph K5", 5, [
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
            (0, 2), (2, 4), (4, 1), (1, 3), (3, 0)
        ], False),
        ("Complete Bipartite K3,3", 6, [
            (0, 3), (0, 4), (0, 5),
            (1, 3), (1, 4), (1, 5),
            (2, 3), (2, 4), (2, 5)
        ], False),
        ("Petersen Graph", 10, [
            (0, 1), (0, 4), (0, 5), (1, 2), (1, 6),
            (2, 3), (2, 7), (3, 4), (3, 8), (4, 9),
            (5, 7), (5, 8), (6, 8), (6, 9), (7, 9)
        ], False),
        ("Hypercube Graph Q3", 8, [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7)
        ], True),
        ("Graph with Self Loops and Parallel Edges", 4, [
            (0, 0), (1, 1), (0, 1), (0, 1), (1, 2), (2, 3), (3, 0)
        ], True),
        ("Disconnected: Two Triangles", 6, [
            (0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)
        ], True),
        ("Disconnected: Triangle + K3,3", 9, [
            (0, 1), (1, 2), (2, 0),
            (3, 6), (3, 7), (3, 8),
            (4, 6), (4, 7), (4, 8),
            (5, 6), (5, 7), (5, 8)
        ], False),
    ]

    passed = 0
    for name, v_count, edges, expected in test_cases:
        tester = LRPlanarityTest(list(range(v_count)), edges)
        res, _, _ = tester.run()
        ok = (res == expected)
        if ok:
            passed += 1
        status = "OK" if ok else "FAIL"
        print(f"[{status}] {name:38} got {res}, expected {expected}")

    print(f"\nPassed {passed}/{len(test_cases)} test cases.")

    # cross-check with networkx if installed
    try:
        import networkx as nx
        import random

        print("\nRunning random graph checks against NetworkX...")
        random.seed(42)
        nx_ok = 0
        total_rnd = 50
        for i in range(total_rnd):
            n = random.randint(3, 16)
            p = random.uniform(0.1, 0.45)
            G = nx.erdos_renyi_graph(n, p, seed=i)
            mapping = {node: idx for idx, node in enumerate(G.nodes())}
            edges = [(mapping[u], mapping[v]) for u, v in G.edges()]
            tester = LRPlanarityTest(list(range(len(G.nodes()))), edges)
            res, _, _ = tester.run()
            if res == nx.is_planar(G):
                nx_ok += 1
        print(f"NetworkX match: {nx_ok}/{total_rnd}")
    except ImportError:
        pass


if __name__ == "__main__":
    run_tests()
