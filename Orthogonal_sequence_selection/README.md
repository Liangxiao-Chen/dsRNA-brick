# Orthogonal (b)KL Sequence Selection

## Big Picture Pipeline

For `N` (default `9`), the workflow has 4 steps:

1. Generate candidate RNA pool
   - Output: `RNA_Pool_Nnt.txt`

2. Filter candidates and build conflict graph
   - Outputs: `FilteredPool_Nnt.txt`, `ConflictGraph_Nnt.txt`

3. Select a large non-conflicting set (independent set)
   - Output: `SelectedPool_from_graph.txt`

4. Generate full complementary pairs and validate format
   - Output: `Orthaganal_RNA_Pool_Nnt.txt`

All files are written into a new folder:

- `Orthaganal_RNA_Pool_Nnt`

(Here `N` is replaced by your `--num-nt` value, e.g. `9`.)

---

## How To Run

```bash
python run_orthogonal_selection.py -N 9
```

Example output folder:

- `Orthaganal_RNA_Pool_9nt/`

Example files inside:

- `RNA_Pool_9nt.txt`
- `FilteredPool_9nt.txt`
- `ConflictGraph_9nt.txt`
- `SelectedPool_from_graph.txt`
- `Orthaganal_RNA_Pool_9nt.txt`

---

## Supported Options

- `-N`, `--num-nt`
  - Meaning: sequence length `N`.
  - Default: `9`.

- `-R`, `--rounds`
  - Meaning: number of randomized attempts in step 3.
  - Default: `1000`.

- `-G`, `--guide`
  - Meaning: guide sequence used in filtering (step 2).
  - Default: `CACGAAGUCAAUAC`.
  - Multiple guides are allowed by repeating `-G`.

- `-O`, `--output`
  - Meaning: filename of step-4 final output.
  - Default: `Orthaganal_RNA_Pool_Nnt.txt`.

- `-U`, `--User`
  - Meaning: show built-in usage instructions.

---

## Command Examples

Use defaults:

```bash
python run_orthogonal_selection.py
```

Set length and rounds:

```bash
python run_orthogonal_selection.py -N 9 -R 2000
```

Use multiple guides:

```bash
python run_orthogonal_selection.py -N 9 -G CACGAAGUCAAUAC -G GGGAAAUUU
```

Set custom final output filename:

```bash
python run_orthogonal_selection.py -N 9 -O My_Final_Pool.txt
```

Show help/instructions:

```bash
python run_orthogonal_selection.py -U
```

---

## Folder Layout

- `run_orthogonal_selection.py`: one-command pipeline runner.
- `functions/generate_pool_nnt.py`: step 1 implementation.
- `functions/build_rna_conflict_graphV2.py`: step 2 implementation.
- `functions/select_from_conflict_graph.py`: step 3 implementation.
- `functions/pool_rna_complement.py`: step 4 implementation.
- `Demo_9nt/`: reference demo outputs.
