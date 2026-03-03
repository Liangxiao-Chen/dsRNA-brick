# dsRNA-brick

This repository provides an end-to-end workflow for dsRNA-brick design:

1. Build an orthogonal KL sequence pool.
2. Build/select lattice structures in 2D or 3D GUI.
3. Generate tile sequences for downstream design workflows.

## Repository Modules

- `Orthogonal_sequence_selection/`
  - KL candidate generation, conflict-graph construction, and independent-set selection.
  - Entry scripts:
    - `generate_pool_nnt.py`
    - `build_rna_conflict_graphV2.py`
    - `select_from_conflict_graph.py`

- `build_dsRNA_bricks_2D/`
  - 2D lattice GUI (`X`, `Y`) with synchronized model/XY map tile selection.
  - Integrated sequence-generation pipeline.

- `build_dsRNA_bricks_3D/`
  - 3D lattice GUI (`X`, `Y`, `Z`) with synchronized 3D/layered-map tile selection.
  - Integrated sequence-generation pipeline.

## Quick Start

### 1. Download the code

Option A (Git):

```bash
git clone https://github.com/Liangxiao-Chen/dsRNA-brick.git dsRNA-brick
cd dsRNA-brick
```

Option B (ZIP):

1. Download ZIP from GitHub.
2. Unzip it.
3. Open terminal in the unzipped `dsRNA-brick` folder.

### 2. Set up Python environment

```bash
conda create -n dsrna python=3.12 -y
conda activate dsrna
pip install "pyside6==6.9.*" numpy pyvista pyvistaqt vtk
```

Install/configure NUPACK separately using official documentation:

- [https://docs.nupack.org/](https://docs.nupack.org/)

### 3. Run orthogonal sequence selection

```bash
cd Orthogonal_sequence_selection
python generate_pool_nnt.py 9
python build_rna_conflict_graphV2.py --num-nt 9 --input RNAPool_9nt.txt
python select_from_conflict_graph.py \
  --seq-file FilteredPool_9nt_out.txt \
  --graph-file ConflictGraph_9nt_edges_out.txt \
  --rounds 1000
```

Main outputs:

- `RNAPool_9nt.txt`
- `FilteredPool_9nt_out.txt`
- `ConflictGraph_9nt_edges_out.txt`
- `SelectedPool_from_graph_out.txt`

### 4. Run the 2D GUI

```bash
cd ../build_dsRNA_bricks_2D
python build_dsRNA_Bricks_2D.py
```

### 5. Run the 3D GUI

```bash
cd ../build_dsRNA_bricks_3D
python build_dsRNA_Bricks_3D.py
```


## Module Documentation

- Orthogonal selection: [`Orthogonal_sequence_selection/README.md`](Orthogonal_sequence_selection/README.md)
- 2D GUI: [`build_dsRNA_bricks_2D/README.md`](build_dsRNA_bricks_2D/README.md)
- 3D GUI: [`build_dsRNA_bricks_3D/README.md`](build_dsRNA_bricks_3D/README.md)

## Repository Structure

```text
dsRNA-brick/
├── Orthogonal_sequence_selection/
├── build_dsRNA_bricks_2D/
├── build_dsRNA_bricks_3D/
├── LICENSE
└── README.md
```

## License

MIT License. See [`LICENSE`](LICENSE).
