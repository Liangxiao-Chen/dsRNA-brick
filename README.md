# dsRNA-brick

An end-to-end workflow for designing **2D dsRNA brick (C-tile) RNA nanostructures**, including:

1) **Orthogonal KL pool selection** (sequence filtering → conflict graph → independent set)  
2) **2D structure design** (GUI-based lattice/shape design and file export)  
3) **NUPACK sequence design** (multi-tile design workflow and result export)

---

## Table of Contents

- [Workflow overview](#workflow-overview)
- [Repository structure](#repository-structure)
- [Quick start (paper example: 9×12)](#quick-start-paper-example-9×12)
- [Module 1 — Orthogonal sequence selection](#module-1--orthogonal-sequence-selection)
- [Module 2 — 2D dsRNA bricks structure design (GUI)](#module-2--2d-dsrna-bricks-structure-design-gui)
- [Module 3 — NUPACK sequence design](#module-3--nupack-sequence-design)
- [NUPACK installation and citation](#nupack-installation-and-citation)
- [Troubleshooting](#troubleshooting)
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

## Quick start (paper example: 9×12)

### 1) Prepare an orthogonal KL pool
Use a final KL pool file generated in:
- `Orthogonal_sequence_selection/`  
Example file: `SelectedPool_from_graph_out.txt`

### 2) Build a 9×12 lattice in the GUI
Run the GUI script inside `build_dsRNA_bricks/` (use the exact filename in that folder):
```bash
python build_dsRNA_Bricks_2D_CTileV3.2.py
```

Then in the GUI:
- Select **KL pool file** = `Orthogonal_sequence_selection/SelectedPool_from_graph_out.txt`
- Set `X = 9`, `Y = 12`, and a `Prefix` (e.g., `Generated_Tiles`)
- Click **Build lattice**
- Click / drag to select tiles defining your target shape
- Click **Generate files** to export:
  - `Generated_Tiles_9x12_108_KL.txt`
  - `Generated_Tiles_9x12_108_NUPACK.txt`
  - `Generated_Tiles_9x12_108.txt`

### 3) Run NUPACK sequence design
Follow the instructions in:
- [`NUPACK_design/README.md`](NUPACK_design/README.md)

---

## Module 1 — Orthogonal sequence selection

Location: `Orthogonal_sequence_selection/`  
Documentation: [`Orthogonal_sequence_selection/README.md`](Orthogonal_sequence_selection/README.md)

Purpose:
- generate candidate KL pools
- remove problematic/self-complementary sequences
- build a conflict graph of overly complementary pairs
- select a large independent set as the final orthogonal pool

---

## Module 2 — 2D dsRNA bricks structure design (GUI)

Location: `build_dsRNA_bricks/`  
Documentation: [`build_dsRNA_bricks/README.md`](build_dsRNA_bricks/README.md)

Purpose:
- interactively design **arbitrary size** and **arbitrary shape** 2D dsRNA brick lattices
- export three downstream files: `*_KL.txt`, `*_NUPACK.txt`, and `*.txt`

Run:
```bash
python build_dsRNA_Bricks_2D_CTileV3.2.py
```

---

## Module 3 — NUPACK sequence design

Location: `NUPACK_design/`  
Documentation: [`NUPACK_design/README.md`](NUPACK_design/README.md)

Purpose:
- run multi-tile NUPACK design for the exported tile definitions
- perform multiple trials per tile and select the best result
- generate consolidated outputs (e.g., `Output_sequence`)

---

## NUPACK installation and citation

Install NUPACK following the official documentation:
- https://docs.nupack.org/

If you use NUPACK in academic work, cite the appropriate NUPACK references as described in the official documentation.

---

## Troubleshooting

- `ModuleNotFoundError: No module named 'nupack'`  
  Ensure NUPACK is installed and your Jupyter kernel matches the environment.

- `len(sequence) != len(structure)`  
  Ensure sequence and dot-paren structure lengths match after removing spaces.

- NUPACK infeasibility (e.g., `No nucleotides found satisfy these constraints`)  
  Reduce soft-constraint weights, or relax discouraged pattern constraints.

---

## License

This project is released under the **MIT License** (see `LICENSE`).

---

### Repo hygiene note (recommended)
If you see `.DS_Store` committed (macOS artifact), remove it and add it to `.gitignore` to prevent future commits:
- Add a line: `.DS_Store`
