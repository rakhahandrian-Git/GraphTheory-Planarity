# Left-Right Planarity Test (Hopcroft-Tarjan Formulation)

## Algorithm Overview

### 1. Algorithm Name & Context
* **Algorithm:** **Left-Right Planarity Test** (Brandes 2009 / Hopcroft & Tarjan 1974)
* **Distinction:**
  * Different from **Kuratowski's theorem** (does not search for $K_5$ / $K_{3,3}$ topological minors or subdivisions).
  * Different from **Narsingh Deo's elementary reduction** (does not simplify via degree-2 vertex contraction or $K_4$ reduction).
  * Directly based on the Left-Right Trémaux tree formulation by de Fraysseix & Rosenstiehl (1982) and simplified by Ulrik Brandes (2009).

### 2. How the Algorithm Works
1. **Euler's Formula Pre-check:**
   For any simple planar graph with $V \ge 3$, the maximum number of edges is $E \le 3V - 6$. If $E > 3V - 6$, the graph is immediately non-planar. Multi-edges and self-loops are cleaned up first.
2. **DFS Tree (Trémaux Tree) Orientation:**
   Run DFS to classify edges into directed **tree edges** and **back-edges** (pointing back to ancestors). Compute vertex heights, return points ($\text{lowpt}$ and $\text{lowpt2}$), and nesting depths.
3. **Nesting Depth Ordering:**
   Sort adjacency lists by nesting depth so paths with lower return points are explored first.
4. **Left-Right Conflict Testing:**
   Back-edges returning to ancestors must be embedded on either the Left or Right side of DFS paths. If two back-edges interlace, they must be assigned to opposite sides.
   The algorithm maintains a stack of conflict pairs (`Interval` objects for left and right return edges). If merging constraints causes a contradiction, the graph is non-planar.
5. **Sign Resolution:**
   If all components pass testing, relative sides are resolved into final Left/Right edge assignments.

### 3. Complexity
* **Time Complexity:** $\mathcal{O}(V + E)$ — linear time, visits each vertex and edge a constant number of times.
* **Space Complexity:** $\mathcal{O}(V + E)$ — storage for adjacency lists, DFS tree heights, lowpoints, and conflict stack.

---

## Prerequisites

* Python 3.8+
* Uses only standard library (`sys`, `argparse`, `collections`). No external packages required to run the planarity test.
* *(Optional)* `networkx` is only used if you run the cross-validation script in `test_planarity.py`.

---

## Instructions to Run

Navigate into this directory:
```bash
cd hopcroft-tarjan
```

### 1. Run with sample inputs file
```bash
python planarity_lr.py input_samples.txt
```

### 2. Run with verbose trace
```bash
python planarity_lr.py input_samples.txt -v
```

### 3. Run built-in demo
```bash
python planarity_lr.py --demo
```

### 4. Run test suite
```bash
python test_planarity.py
```

---

## Result of Sample Run

Output from running `python planarity_lr.py input_samples.txt`:

```text
Found 10 graph(s) to evaluate.
--------------------------------------------------

>> Graph: 'K4' (V=4, E=6)
   Verdict: [ PLANAR ]
   Details: Graph is planar
   LR-Partition (6 edges):
     * 0 -> 3: Tree Edge (Right)
     * 2 -> 1: Back Edge (Right)
     * 1 -> 0: Tree Edge (Right)
     * 3 -> 1: Back Edge (Right)
     * 2 -> 0: Back Edge (Right)
     * 3 -> 2: Tree Edge (Right)

>> Graph: 'K5' (V=5, E=10)
   Verdict: [ NON-PLANAR ]
   Details: Euler formula check failed: E (10) > 3V - 6 (9)

>> Graph: 'K3_3' (V=6, E=9)
   Verdict: [ NON-PLANAR ]
   Details: Conflict detected: back-edges cannot be partitioned without crossing

>> Graph: 'Cube_Q3' (V=8, E=12)
   Verdict: [ PLANAR ]
   Details: Graph is planar
   LR-Partition (12 edges):
     * 7 -> 3: Back Edge (Right)
     * 4 -> 7: Tree Edge (Right)
     * 6 -> 5: Tree Edge (Right)
     * 5 -> 1: Back Edge (Right)
     * 3 -> 2: Tree Edge (Right)
     * 1 -> 0: Tree Edge (Right)
     * 2 -> 1: Back Edge (Right)
     * 4 -> 0: Back Edge (Right)
     * 2 -> 6: Tree Edge (Right)
     * 7 -> 6: Back Edge (Right)
     * 5 -> 4: Tree Edge (Right)
     * 0 -> 3: Tree Edge (Right)

>> Graph: 'Petersen' (V=10, E=15)
   Verdict: [ NON-PLANAR ]
   Details: Conflict detected: back-edges cannot be partitioned without crossing

>> Graph: 'Pentagonal_Prism' (V=10, E=15)
   Verdict: [ PLANAR ]
   Details: Graph is planar
   LR-Partition (15 edges):
     * 7 -> 6: Tree Edge (Right)
     * 2 -> 7: Tree Edge (Right)
     * 3 -> 2: Tree Edge (Right)
     * 8 -> 7: Back Edge (Right)
     * 9 -> 4: Back Edge (Right)
     * 5 -> 9: Tree Edge (Right)
     * 2 -> 1: Back Edge (Right)
     * 4 -> 3: Tree Edge (Right)
     * 1 -> 0: Tree Edge (Right)
     * 6 -> 1: Back Edge (Right)
     * 5 -> 0: Back Edge (Right)
     * 9 -> 8: Tree Edge (Right)
     * 8 -> 3: Back Edge (Right)
     * 0 -> 4: Tree Edge (Right)
     * 6 -> 5: Tree Edge (Right)

>> Graph: 'Lecture_Slide4' (V=7, E=11)
   Verdict: [ PLANAR ]
   Details: Graph is planar
   LR-Partition (11 edges):
     * 1 -> 0: Tree Edge (Right)
     * 3 -> 2: Tree Edge (Right)
     * 4 -> 3: Tree Edge (Right)
     * 6 -> 1: Back Edge (Right)
     * 0 -> 6: Tree Edge (Right)
     * 6 -> 5: Tree Edge (Right)
     * 2 -> 1: Back Edge (Right)
     * 5 -> 4: Tree Edge (Right)
     * 2 -> 5: Back Edge (Right)
     * 3 -> 5: Back Edge (Right)
     * 2 -> 6: Back Edge (Right)

>> Graph: 'Lecture_Slide14' (V=8, E=12)
   Verdict: [ PLANAR ]
   Details: Graph is planar
   LR-Partition (12 edges):
     * 7 -> 6: Tree Edge (Right)
     * 4 -> 0: Tree Edge (Right)
     * 9 -> 6: Back Edge (Right)
     * 6 -> 8: Tree Edge (Right)
     * 0 -> 5: Back Edge (Right)
     * 5 -> 4: Tree Edge (Right)
     * 0 -> 2: Tree Edge (Right)
     * 9 -> 7: Back Edge (Right)
     * 2 -> 4: Back Edge (Right)
     * 8 -> 7: Back Edge (Right)
     * 2 -> 5: Back Edge (Right)
     * 8 -> 9: Tree Edge (Right)

>> Graph: 'Wheel_W6' (V=6, E=10)
   Verdict: [ PLANAR ]
   Details: Graph is planar

>> Graph: 'Disconnected_K3_3' (V=9, E=12)
   Verdict: [ NON-PLANAR ]
   Details: Conflict detected: back-edges cannot be partitioned without crossing
```

---

## AI Tools Usage Disclosure

Per slide 19 & 20 guidelines:
* **Tool:** Google Gemini
* **Scope of usage:**
  * Referencing Ulrik Brandes' 2009 Left-Right planarity formulation and data structures.
  * Writing test cases and cross-validation script.
  * Assisting with Markdown formatting.
* **Prompt history summary:**
  * Provided lecture PDF and requested implementation of 1 planarity testing algorithm (excluding Kuratowski and Deo's reduction) with test cases and README.
  * Adjusted code style and comments to be concise, natural, and readable.

---

## References
1. Hopcroft, J., & Tarjan, R. (1974). *Efficient Planarity Testing.* JACM.
2. de Fraysseix, H., & Rosenstiehl, P. (1982). *A characterization of planar graphs by Trémaux trees.* Combinatorica.
3. Brandes, U. (2009). *The Left-Right Planarity Test.*
4. Deo, N. (1974). *Graph Theory with Applications to Engineering and Computer Science.*
5. AI lol
