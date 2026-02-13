# dsRNA-brick

An end-to-end workflow for designing **2D dsRNA brick (C-tile) RNA nanostructures**, including:

1) **Orthogonal KL pool selection** (sequence filtering → conflict graph → independent set)  
2) **2D structure design** (GUI-based lattice/shape design and file export)  
3) **NUPACK sequence design** (multi-tile design workflow and result export)

---

## Table of Contents

- [Workflow overview](#workflow-overview)
- [Repository structure](#repository-structure)
- [License](#license)

---

## Workflow overview

### Step 1 — Orthogonal sequence selection
Generate a candidate KL pool, remove self-/guide-complementary sequences, build a **conflict graph**, and select a **large independent set** as the final orthogonal KL pool.

**Typical output:**
- `SelectedPool_from_graph_out.txt` (final orthogonal KL pool)

### Step 2 — 2D structure design (GUI)
Use a GUI to define the lattice size (**X cols**, **Y rows**) and interactively select an arbitrary 2D shape (tiles). Export design files for downstream sequence design and documentation.

**Typical outputs:**
- `*_KL.txt` — KL pairing relationships in the generated structure  
- `*_NUPACK.txt` — input text for NUPACK sequence design  
- `*.txt` — per-tile information (tile geometry and secondary structure)

### Step 3 — NUPACK sequence design
Run NUPACK-based multi-tile sequence design using the exported structure/tile definitions. The workflow runs multiple trials per tile, selects the best result, and writes a consolidated output (e.g., `Output_sequence`).

---

## Repository structure

- `Orthogonal_sequence_selection/`  
  Scripts and outputs for KL pool generation, conflict-graph construction, and orthogonal pool selection.  
  See: [`Orthogonal_sequence_selection/README.md`](Orthogonal_sequence_selection/README.md)

- `build_dsRNA_bricks/`  
  GUI tool for designing 2D dsRNA brick lattices and exporting structure files for downstream design.  
  See: [`build_dsRNA_bricks/README.md`](build_dsRNA_bricks/README.md)

- `NUPACK_design/`  
  NUPACK multi-tile design workflow (notebooks/scripts + documentation).  
  See: [`NUPACK_design/README.md`](NUPACK_design/README.md)

---

## License

This project is released under the **MIT License** (see `LICENSE`).

---

### Repo hygiene note (recommended)
If you see `.DS_Store` committed (macOS artifact), remove it and add it to `.gitignore` to prevent future commits:
- Add a line: `.DS_Store`
