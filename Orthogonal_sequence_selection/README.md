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

Goal:

- Enumerate “good” `N`-nt RNA candidates with controlled GC/AU composition and sequence-quality constraints.

Core idea:

- Build candidates from strong/weak patterns:
  - `S` (strong) = `G/C`
  - `W` (weak) = `A/U`
- For `N=9`, the generator uses exactly 4 `S` and 5 `W`.

Major filters used in this step:

- Avoid problematic S/W runs (e.g., long GC-only stretches).
- Keep all four nucleotides represented (`A`, `U`, `C`, `G`).
- Hairpin-risk exclusion near ends for `N >= 9`.
- Exclude internal `UUU` / mirrored `AAA` patterns.
- Additional end-quality checks and degenerate-run checks (`K/M/R/Y`).

Final output:

- `RNA_Pool_Nnt.txt`

---

### Step 2: Filter pool and build conflict graph

Goal:

- Remove sequences that are unsuitable for orthogonal KL use, then build a conflict graph capturing pairwise incompatibility.

Pre-filtering:

- Remove self-complementary (palindromic) sequences.
- Remove sequences that are too complementary to guide sequence(s) (`-G`), by scanning all contiguous `N`-mers of each guide.

Conflict graph construction:

- One node per remaining sequence.
- Undirected edge = two sequences are too complementary under the script criteria.
- Pairing model includes:
  - Watson-Crick pairs (`A-U`, `U-A`, `G-C`, `C-G`)
  - GU wobble (`G-U`, `U-G`)
- Checks both direct and reverse-complement orientations.
- Includes bulged-pair checks for longer sequences (`N >= 9`).

Final outputs:

- `FilteredPool_Nnt.txt`
- `ConflictGraph_Nnt.txt`

---

### Step 3: Select a large independent set

Goal:

- Find a large subset of sequences with no conflicts between any pair.

Method:

- Randomized greedy independent-set search for `-R` rounds.
- Keep the largest set found across rounds.
- Run an additional optimization pass:
  - remove one seed element,
  - regrow greedily,
  - keep improvements.

Interpretation:

- The result is a large set of mutually non-conflicting candidates from the graph model.

Final output:

- `SelectedPool_from_graph.txt`

---

### Step 4: Generate full complementary pool and validate format

Goal:

- Convert selected sequences into final KL-pool lines (`sequence` + full complement) and ensure downstream-compatible format.

What this step does:

- Reads selected entries from step 3.
- Generates full complementary sequence for each selected sequence.
- Writes tab-separated pair lines:
  - `sequence<TAB>complement_sequence`
- Validates output format for downstream tile-design pipelines (one pair per non-empty line).

Final output:

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
