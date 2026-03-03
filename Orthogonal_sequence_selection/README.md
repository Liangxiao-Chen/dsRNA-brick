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

## Detailed Step Introduction

### Step 1: Generate candidate RNA pool

This step enumerates candidate `N`-nt RNA sequences with the intended base-composition and local-structure constraints.

Main logic:

- Uses strong/weak patterning (`S`: `G/C`, `W`: `A/U`) to control composition.
- Applies sequence-quality filters (e.g., avoid problematic end patterns and long degenerate runs).
- Removes obvious hairpin-prone/internal unstable patterns.

Final output of step 1:

- `RNA_Pool_Nnt.txt`

### Step 2: Filter pool and build conflict graph

This step removes problematic candidates and then builds an undirected conflict graph.

Main logic:

- Removes self-complementary (palindromic) sequences.
- Removes sequences that strongly conflict with guide sequence(s) (`-G`).
- Compares remaining sequences pairwise under Watson-Crick + GU wobble rules.
- Adds a graph edge for pairs considered too complementary.

Final outputs of step 2:

- `FilteredPool_Nnt.txt`
- `ConflictGraph_Nnt.txt`

### Step 3: Select a large independent set

This step selects a large non-conflicting subset from the conflict graph.

Main logic:

- Runs randomized greedy independent-set search for `-R` rounds.
- Keeps the best result found across rounds.
- Applies an additional optimization pass to try improving the selected set.

Final output of step 3:

- `SelectedPool_from_graph.txt`

### Step 4: Generate full complementary pool and validate format

This step converts selected sequences into complementary pairs and writes the final pool file for downstream tile-design tools.

Main logic:

- Reads selected entries from step 3.
- Writes `sequence<TAB>complement_sequence` pairs.
- Ensures final file format is compatible with downstream KL-pool input requirements.

Final output of step 4:

- `Orthaganal_RNA_Pool_Nnt.txt`

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
