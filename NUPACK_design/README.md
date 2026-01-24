# NUPACK RNA Tile Sequence Design (Jupyter / Multi-tile)

This repository provides a Jupyter notebook workflow for **large-scale RNA sequence design** using **NUPACK**.  
It reads an input file containing many `TILE` blocks (sequence + dot-paren secondary structure), runs NUPACK design trials per tile, selects the best result, and writes a consolidated output file.

## Installation (NUPACK)

Install NUPACK by following the official documentation:

https://docs.nupack.org/

After installation, ensure the Jupyter kernel you use can import `nupack`.

## What the notebook does

For each tile in the input file, the notebook will:

- Parse the `*******TILEi*******` header, sequence line, and structure line
- Remove spaces for computation (spaces are treated as human-readable separators)
- Validate `len(sequence) == len(structure)` after cleaning
- Run **NUPACK `tube_design`** with a configurable number of trials (default: `TRIALS = 3`)
- Choose the **best trial** by **lowest ensemble defect**
- Restore the original spacing/segmentation in the reported sequences for readability
- Append results for all tiles into one output file named **`Output_sequence`**

The design uses **soft constraints** (pattern penalties) to discourage undesired motifs while reducing the risk of infeasible designs.

## Input file format

Plain text file with an optional header line, followed by repeated tile blocks:

```text
(optional header line)
*******TILE0*******
<sequence line, may contain spaces>
<structure line, dot-paren, may contain spaces>

*******TILE1*******
<sequence line>
<structure line>
...
```

Notes:

- The notebook **removes spaces** before calling NUPACK.
- The secondary structure must match the sequence length **after spaces are removed**.
- Any `T` in the sequence is automatically converted to `U`.

## Output file format

The output file is named:

```text
Output_sequence
```

It contains one block per tile:

```text
(optional header line)

*******TILE0*******

INPUT_SEQUENCE	<sequence with original spacing restored>

INPUT_STRUCTURE	
<structure with original spacing restored>

OUTPUT_SEQUENCE	<best designed sequence with original spacing restored>
```

The output is intentionally minimal and does **not** repeat the raw input tile lines beyond the tile header.

## How to run

1. Open the notebook (`.ipynb`) in Jupyter Lab / Jupyter Notebook.
2. Place the input file in the same directory as the notebook (or update the path).
3. In the config cell, set:

```python
INPUT_FILE = "Generated_Tiles_9x12_108.txt"  # or your filename
TRIALS = 3
```

4. Run all cells.
5. The results will be written to `Output_sequence`.

## Soft-constraint weight (pattern penalties)

Soft constraints use a `weight` parameter:

- Higher `weight` → stronger preference to avoid the specified patterns, but higher risk of infeasibility
- Lower `weight` → easier feasibility, but more pattern violations may appear

If you encounter an error like:

- `No nucleotides found satisfy these constraints`

reduce the soft-constraint `weight` (e.g., `0.25 → 0.1 → 0.05`), or reduce the number of discouraged patterns.

## Troubleshooting

### `ModuleNotFoundError: No module named 'nupack'`
- Confirm NUPACK is installed (per the official docs).
- Confirm your Jupyter kernel points to the same environment where NUPACK is installed.

### Length mismatch errors
If the notebook reports `len(sequence) != len(structure)`:
- Ensure the structure line includes exactly one character per nucleotide (dot-paren), ignoring spaces.
- Verify you did not accidentally remove or add nucleotides in the sequence line.

## Citation

If you use NUPACK in academic work, cite the appropriate NUPACK references as described in the official documentation:

https://docs.nupack.org/
