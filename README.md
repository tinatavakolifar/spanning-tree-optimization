# Spanning Tree Optimization
## Overview

This project demonstrates the computation and visualization of several spanning tree optimizations from a weighted undirected graph loaded from an adjacency matrix file named `adj_matrix.csv`.

## Spanning Tree Types Covered

- **MST**: Minimum Spanning Tree (minimizes total weight)  
- **MBST**: Minimum Bottleneck Spanning Tree (minimizes the largest edge)  
- **MMST**: Minimum Median Spanning Tree (minimizes the median edge weight)  
- **MVST**: Minimum Variance Spanning Tree (minimizes weight variance)  

All results are computed, visualized, and explained in a Jupyter notebook.

---

## Repository Contents

- `Spanning_Tree_OP.py` — Main Python script, defines all functions, classes, and exports results to PDF  
- `Spanning_Tree_OP.ipynb` — Fully documented Jupyter Notebook with explanations and step-by-step visualizations  
- `adj_matrix.csv` — Input adjacency matrix (symmetric, weighted)   
- `intro_to_MBST_MMST_MVST.pdf` - Theoretical overview of MST, MBST, MMST and MVST

---

## Installation

1. **Clone the repository:**  
```bash
git clone https://github.com/tinatavakolifar/Spanning_Tree_OP.git
cd spanning-tree-optimization
```

2. **Requirements**

Python ≥ 3.9

Libraries: numpy ≥ 1.25, matplotlib ≥ 3.8, networkx ≥ 3.1

install dependencies with:
```bash
pip install -r requirements.txt
```

3. **Run the Jupyter notebook:**  

```bash
jupyter notebook Spanning_Tree_OP.ipynb
```
or run the Python script:
```bash
python SameMSTfile_but_in.py
```

---

## Visual Outputs

The project generates visualizations for:

1. Original Graph
2. Minimum Spanning Tree (MST)
3. Minimum Bottleneck Spanning Tree (MBST)
4. Minimum Median Spanning Tree (MMST)
5. Minimum Variance Spanning Tree (MVST)

All figures are also saved in `graph_report.pdf`.

---

## Algorithms Used

- Kruskals Algorithm — finds MST efficiently
- Disjoint Set Union (DSU) — detects and prevents cycles
- Edge Replacement Strategy — explores near-MST variants to find MBST, MMST, MVST
- NetworkX + Matplotlib — visualization

---

## Purpose / Learning Outcomes

This project explores how different definitions of spanning trees affect graph structure and the trade-offs between total weight, edge uniformity, and bottleneck reduction. It also provides practical understanding of Kruskals algorithm and its variations.

---

## Credits

Author: Tina Tavakolifar  
