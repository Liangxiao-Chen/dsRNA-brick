# Orthogonal (b)KL sequence selection

> Source: Evernote export (cleaned for GitHub Markdown)

## Big picture: what the three scripts do

For **N = 9** (our standard KL length), the pipeline is:

1. **`generate_pool_nnt.py`**  
   Generate a *raw* pool of “good” 9-nt sequences with nice local properties (GC%, homopolymer limit, avoid specific motifs, etc.).  
   Output: a text file of candidate sequences.

2. **`orthogonal_set_from_pool.py`**  
   From the raw pool, greedily / iteratively pick a subset that satisfies global orthogonality constraints (pairwise similarity, mismatch/shift checks, etc.).  
   Output: an “orthogonal set” list.

3. **`pick_kl_to_kd.py`**  
   Given the orthogonal KL list, generate the corresponding KD set (or mapping) with any additional constraints needed for your brick design.  
   Output: final mapping / sequences for downstream design.

---

## Script 1: `generate_pool_nnt.py`

### Summary
Generate a pool of candidate **N-nt** sequences (default **N=9**) filtered by local sequence rules.

### Typical usage
```bash
python generate_pool_nnt.py   --n 9   --gc_min 0.33 --gc_max 0.67   --max_homopolymer 3   --out pool_9nt.txt
```

### Key filters (typical)
- GC content range
- Maximum homopolymer length (e.g., `AAAA` not allowed)
- Avoid specific motifs / restriction sites
- Optional: avoid hairpin-prone patterns (if implemented)

### Output
A plain text file, one sequence per line, e.g.
```text
ACGTTGCAA
GCTACATGG
...
```

---

## Script 2: `orthogonal_set_from_pool.py`

### Summary
Select an orthogonal subset from the candidate pool with global constraints.

### Typical usage
```bash
python orthogonal_set_from_pool.py   --pool pool_9nt.txt   --set_size 48   --max_pair_identity 5   --out orthogonal_kl_48.txt
```

### Notes on orthogonality checks (common examples)
- Pairwise identity threshold
- Shifted alignment checks (to avoid off-register binding)
- Optional: reverse-complement checks
- Optional: energy-based screening (if coupled to a model)

### Output
A list of selected KL sequences, e.g.
```text
KL01  ACGTTGCAA
KL02  GCTACATGG
...
```

---

## Script 3: `pick_kl_to_kd.py`

### Summary
Map KL sequences to KD sequences (and/or generate KD set) for the dsRNA-brick design workflow.

### Typical usage
```bash
python pick_kl_to_kd.py   --kl orthogonal_kl_48.txt   --out kl_kd_map.csv
```

### Output
Typically a CSV/TSV mapping table, e.g.
```text
KL_id,KL_seq,KD_id,KD_seq
KL01,ACGTTGCAA,KD01,TTGCAACGT
...
```

---

## Practical tips / troubleshooting
- If the pool is too small, relax filters in `generate_pool_nnt.py` (e.g., widen GC range, increase homopolymer limit).
- If orthogonal selection fails or is slow, reduce target set size or loosen pairwise thresholds in `orthogonal_set_from_pool.py`.
- Keep outputs small in the repo; if you have large pools/sets, store them as **GitHub Releases** or external links, and commit only minimal examples.

---

## Suggested repo organization
- `src/` core algorithm code
- `gui/` GUI app for interactive design
- `examples/` example inputs/outputs
- `docs/` detailed documentation (recommended place for extended notes)
