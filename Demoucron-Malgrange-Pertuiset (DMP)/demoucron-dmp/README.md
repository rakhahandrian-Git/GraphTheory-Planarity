# Demoucron–Malgrange–Pertuiset (DMP) Planarity Test

## Algorithm Overview

### 1. Algorithm Name & Context
* **Algorithm:** **Demoucron–Malgrange–Pertuiset (DMP)** Planarity Test (Demoucron, Malgrange, & Pertuiset, 1964).
* **Classification:** Incremental cycle/path addition (bridge-based) planarity testing.
* **Distinction:**
  * **Different from Kuratowski's Theorem:** It does **not** search for subdivisions of Kuratowski's subgraphs ($K_5$ or $K_{3,3}$).
  * **Different from Narsingh Deo's Elementary Reduction:** It does **not** contract degree-2 vertices or reduce graphs to $K_4$.
  * Instead, it starts with an embedded cycle and iteratively adds paths from un-embedded "bridges" (fragments) into admissible faces of the current planar embedding.

---

### 2. How the Algorithm Works

1. **Preprocessing & Euler's Formula:**
   * Removes self-loops and parallel edges (neither impacts planarity).
   * For any simple planar graph with $V \ge 3$, the inequality $E \le 3V - 6$ must hold. If $E > 3V - 6$, the graph is immediately declared **non-planar**.
   * Decomposes the graph into connected components and **biconnected components (blocks)** using Tarjan's bridge/cutpoint DFS algorithm. A graph is planar if and only if every block is planar.

2. **Initial Face Embedding:**
   * In each 2-connected block ($V \ge 3$), DFS finds an initial simple cycle $C$.
   * $C$ divides the plane/sphere into exactly two initial faces: an interior face (clockwise traversal of $C$) and an exterior face (counter-clockwise traversal of $C$).
   * $H$ stores the current embedded subgraph (initially $H = C$).

3. **Bridge (Fragment) Decomposition:**
   * Relative to $H$, the remaining un-embedded edges ($E \setminus E(H)$) are partitioned into **bridges**:
     * **Chords:** Single edges whose both endpoints already lie on $H$.
     * **Component Bridges:** Connected components of vertices not yet in $H$, along with all their incident edges connecting to $H$.
   * The vertices of $H$ connected to a bridge are its **attachments**.

4. **Admissible Faces & Planarity Check:**
   * For every bridge $B$, determine which faces of the current embedding contain **all** attachment vertices of $B$ along their boundary cycle ($F(B)$).
   * **Planarity Criterion:** If any bridge $B$ has $|F(B)| = 0$, that bridge cannot be embedded into any face without crossing existing edges $\implies$ the graph is **NON-PLANAR**.

5. **Forced Selection & Path Addition ($\alpha$-Path):**
   * If there exists a bridge with $|F(B)| = 1$, it is a **forced choice** and must be embedded immediately into its sole admissible face.
   * Otherwise, select any bridge with minimum $|F(B)|$.
   * Inside the selected bridge, find an $\alpha$-path $P$ connecting two distinct attachment vertices $u$ and $v$ through the bridge.
   * Embed $P$ into the target face $F$, which splits $F$ into two new simple cyclic faces.
   * Add $P$ to $H$ and update the set of faces.

6. **Completion:**
   * Repeat until all edges are embedded. If all biconnected blocks succeed, the graph is **PLANAR**, and the exact list of planar face boundaries is returned.

---

### 3. Complexity
* **Time Complexity:** $\mathcal{O}(V^2)$ — computing bridges and finding admissible faces takes linear time per step; at most $\mathcal{O}(V)$ paths are added.
* **Space Complexity:** $\mathcal{O}(V + E)$ — linear memory for adjacency lists, bridge tracking, and boundary cycles of planar faces.

---

## Prerequisites

* Python 3.8+ (tested up to Python 3.14)
* Standard library only (`sys`, `argparse`, `collections`). No external packages are required to run the planarity test.
* *(Optional)* `networkx` is used only if you wish to run the randomized cross-validation suite in `test_planarity.py`.

---

## Instructions to Run

Navigate into this directory:
```bash
cd demoucron-dmp
```

### 1. Run with Sample Inputs File
```bash
python planarity_dmp.py input_samples.txt
```

### 2. Run with Verbose Step-by-Step Trace
```bash
python planarity_dmp.py input_samples.txt -v
```

### 3. Run Built-in Demo
```bash
python planarity_dmp.py --demo
```

### 4. Run Test Suite
```bash
python test_planarity.py
```

---

## Result of Sample Run

Output from executing `python planarity_dmp.py input_samples.txt`:

```text
Found 10 graph(s) to evaluate.
--------------------------------------------------

>> Graph: 'K4' (V=4, E=6)
   Verdict: [ PLANAR ]
   Details: Graph is planar (all biconnected blocks successfully embedded)
   Planar Faces (4 faces):
     * Face 1: 0 -> 1 -> 3 -> 0
     * Face 2: 1 -> 0 -> 2 -> 1
     * Face 3: 3 -> 1 -> 2 -> 3
     * Face 4: 2 -> 0 -> 3 -> 2

>> Graph: 'K5' (V=5, E=10)
   Verdict: [ NON-PLANAR ]
   Details: Euler formula check failed: E (10) > 3V - 6 (9)

>> Graph: 'K3_3' (V=6, E=9)
   Verdict: [ NON-PLANAR ]
   Details: Block 1/1 non-planar: Bridge with attachments ['1', '0', '2'] has 0 admissible faces (cannot be drawn without crossing)

>> Graph: 'Cube_Q3' (V=8, E=12)
   Verdict: [ PLANAR ]
   Details: Graph is planar (all biconnected blocks successfully embedded)
   Planar Faces (6 faces):
     * Face 1: 5 -> 4 -> 0 -> 1 -> 5
     * Face 2: 4 -> 5 -> 6 -> 7 -> 4
     * Face 3: 6 -> 5 -> 1 -> 2 -> 6
     * Face 4: 0 -> 4 -> 7 -> 3 -> 0
     * Face 5: 3 -> 7 -> 6 -> 2 -> 3
     * Face 6: 2 -> 1 -> 0 -> 3 -> 2

>> Graph: 'Petersen' (V=10, E=15)
   Verdict: [ NON-PLANAR ]
   Details: Block 1/1 non-planar: Bridge with attachments ['6', '1', '3', '5', '7'] has 0 admissible faces (cannot be drawn without crossing)

>> Graph: 'Pentagonal_Prism' (V=10, E=15)
   Verdict: [ PLANAR ]
   Details: Graph is planar (all biconnected blocks successfully embedded)
   Planar Faces (7 faces):
     * Face 1: 6 -> 7 -> 2 -> 1 -> 6
     * Face 2: 6 -> 1 -> 0 -> 5 -> 6
     * Face 3: 2 -> 7 -> 8 -> 3 -> 2
     * Face 4: 8 -> 7 -> 6 -> 5 -> 9 -> 8
     * Face 5: 0 -> 1 -> 2 -> 3 -> 4 -> 0
     * Face 6: 4 -> 3 -> 8 -> 9 -> 4
     * Face 7: 9 -> 5 -> 0 -> 4 -> 9

>> Graph: 'Lecture_Slide4' (V=7, E=11)
   Verdict: [ PLANAR ]
   Details: Graph is planar (all biconnected blocks successfully embedded)
   Planar Faces (6 faces):
     * Face 1: 6 -> 1 -> 0 -> 6
     * Face 2: 1 -> 6 -> 2 -> 1
     * Face 3: 2 -> 6 -> 5 -> 2
     * Face 4: 2 -> 5 -> 3 -> 2
     * Face 5: 3 -> 5 -> 4 -> 3
     * Face 6: 5 -> 6 -> 0 -> 1 -> 2 -> 3 -> 4 -> 5

>> Graph: 'Lecture_Slide14' (V=8, E=12)
   Verdict: [ PLANAR ]
   Details: Graph is planar (all biconnected blocks successfully embedded)
   Planar Faces (8 faces):
     * Face 1: 6 -> 8 -> 9 -> 6
     * Face 2: 8 -> 6 -> 7 -> 8
     * Face 3: 9 -> 8 -> 7 -> 9
     * Face 4: 7 -> 6 -> 9 -> 7
     * Face 5: 4 -> 0 -> 5 -> 4
     * Face 6: 0 -> 4 -> 2 -> 0
     * Face 7: 5 -> 0 -> 2 -> 5
     * Face 8: 2 -> 4 -> 5 -> 2

>> Graph: 'Wheel_W6' (V=6, E=10)
   Verdict: [ PLANAR ]
   Details: Graph is planar (all biconnected blocks successfully embedded)
   Planar Faces (6 faces):
     * Face 1: 0 -> 4 -> 5 -> 0
     * Face 2: 5 -> 1 -> 0 -> 5
     * Face 3: 5 -> 4 -> 3 -> 5
     * Face 4: 3 -> 4 -> 0 -> 1 -> 2 -> 3
     * Face 5: 5 -> 3 -> 2 -> 5
     * Face 6: 2 -> 1 -> 5 -> 2

>> Graph: 'Disconnected_K3_3' (V=9, E=12)
   Verdict: [ NON-PLANAR ]
   Details: Block 1/2 non-planar: Bridge with attachments ['4', '3', '5'] has 0 admissible faces (cannot be drawn without crossing)
```

---

## AI Tools Usage Disclosure

In compliance with the assignment rules, AI assistance (Google Antigravity / Gemini) was utilized during the development of this implementation:
* **Assistance Provided:** Structuring the biconnected component decomposition and boundary cycle slicing for planar face subdivision in the Demoucron-Malgrange-Pertuiset algorithm, as well as formulating automated verification against `networkx.is_planar`.
* **Prompt History:**
  1. *"what is the homework what am i supposed to do"*
  2. *"https://github.com/rakhahandrian-Git/GraphTheory-Planarity can you open this?"*
  3. *"yes choose algorithm that is the easy one"*
