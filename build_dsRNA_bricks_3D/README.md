# dsRNA Bricks 3D GUI (V5_15)

## Overview
This project provides a Qt + PyVista GUI to:

1. Build a 3D dsRNA C-tile lattice (`L` and `W` tile families).
2. Select/deselect tiles in 3D and 2D views.
3. Generate tile sequence + secondary structure blocks for selected tiles.
4. Run NUPACK sequence design on generated tiles, with a progress dialog.

Current version: `V5_15`.

## Project Layout
- `build_dsRNA_Bricks_3D.py`: top-level launcher.
- `function/app.py`: main GUI window and workflow wiring.
- `function/lattice_builder.py`: lattice/tile placement rules.
- `function/map2d.py`: 2D XY layer map widgets and selection.
- `function/view3d.py`: embedded PyVista 3D renderer.
- `function/c_tiles.py`: tile geometry definitions.
- `function/rna_tile_generator.py`: Type I / Type II RNA scaffold generation.
- `function/nupack_runner.py`: NUPACK run pipeline for generated tiles.
- `SelectedPool_from_graph_out_complement.txt`: example KL pool input.
- `NUPACK_sequence_design.ipynb`: original notebook reference.
- `previous/`: archived older versions and historical assets.

## Requirements
Tested stack:
- Python `3.12`
- `PySide6` (`6.9.*` recommended)
- `numpy`
- `pyvista`
- `pyvistaqt`
- `vtk`
- `nupack` (for NUPACK design step)

Example install (inside a conda env):

```bash
conda create -n dsrna python=3.12 -y
conda activate dsrna
pip install "pyside6==6.9.*" numpy pyvista pyvistaqt vtk nupack
```

## Run
From this folder:

```bash
python build_dsRNA_Bricks_3D.py
```

Alternative:

```bash
python -m function.build_dsRNA_Bricks_3D
```

## GUI Workflow
1. Choose KL pool file (`Browse...`).
2. Set `X`, `Y`, `Z` and optional prefix.
3. Click `Build lattice`.
4. Select tiles (3D click or 2D layer maps).
5. Click `Generate sequence` (this now exports 2D SVG and runs NUPACK automatically; progress dialog is shown).

## Sequence + Export + NUPACK Behavior
- Generation uses **selected tiles only**.
- Validation checks include:
  - KL pool availability.
  - Structure closure constraints.
  - Selected-tile connectivity.
- `Generate sequence` now creates three files in the current folder:
  - `<prefix>_tile_structure.txt`
  - `<prefix>_2D_map.svg`
  - `<prefix>_sequence.txt` (NUPACK output)
- Progress dialog remains visible during NUPACK design.

## Notes
- `.gitignore` excludes generated outputs, cache files, and `.DS_Store`.
- `previous/` is intentionally kept for historical reference.
