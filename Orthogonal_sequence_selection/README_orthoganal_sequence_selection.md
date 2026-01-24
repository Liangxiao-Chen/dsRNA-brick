# Orthogonal (b)KL sequence selection

## Big picture: what the three scripts do

For N = 9 (our standard KL length), the pipeline is:

1. `generate_pool_nnt.py` ([Download link](https://www.dropbox.com/scl/fi/0xktm97phin3rtkx8ouf0/generate_pool_nnt.py?rlkey=8skbtq0iathk69ybi80q0tybd&dl=0))
   Generate a **raw** pool of “good” 9‑nt sequences with nice local properties
   → `RNAPool_9nt.txt`
2. `build_rna_conflict_graphV2.py` ([Download link](Orthogonal_sequence_selection/build_rna_conflict_graphV2.py))
   Remove self‑complementary / guide‑complementary sequences and build a **conflict graph** that records which sequences are too complementary to each other
   → `FilteredPool_9nt_out.txt` + `ConflictGraph_9nt_edges_out.txt`
3. `select_from_conflict_graph.py` ([Download link](https://www.dropbox.com/scl/fi/meeqmmzmhzl5t8m1jdzp5/select_from_conflict_graph.py?rlkey=onfc1v6x1j3q5m72xn1b8w5q4&dl=0))
   On that graph, find a **large independent set** (a big subset of sequences that **do not** conflict with each other)
   → final orthogonal set, e.g. `SelectedPool_from_graph_out.txt`

That final set can then be fed into the tile‑building scripts. 

Below I’ll walk through each script (how to run it and what it’s actually doing), then summarize the 9‑nt design criteria and the graph‑theory picture.



---

## Step 1 – `generate_pool_nnt.py`

**Goal:** enumerate “nice” N‑nt KL candidates (for us N = 9) with controlled GC content, no obvious hairpins, and reasonable end patterns.

### How to run it (9‑nt example)

From the project directory:


```
python generate_pool_nnt.py 9
```
- **Positional argument**: `N` (here `9`)
- **Optional**: `--output <filename>`

If you don’t give `--output`, it will write to: `RNAPool_9nt.txt` (one sequence per line, all uppercase RNA).

### What it **means** by “4S and (N–4) W”

- It works at two levels:
1. **S/W pattern** level: "strong" vs "weak"
- **S = strong** = G or C
- **W = weak**  = A or U

3. **Base assignment** level: choose actual A/G/U/C consistent with that S/W pattern.

For N = 9 the script enforces **exactly 4 S and 5 W**, so **every 9‑mer has 4 G/C and 5 A/U**.

#### S/W pattern generation (function `gen_pattern_sw`)

For each 9‑nt pattern, it:

- Picks positions of the 4 S’s (indices Sa, Sb, Sc, Sd) and fills the rest with W.
- Rejects patterns that:
- have **4 consecutive S** (avoids long GC runs),
- have a **long run of W at the 5′ end** (avoids very A/U‑rich 5′ segments, following the original C++ logic).

Result: a list of strings like `"swwswwsww"` that tell you where G/C vs A/U are allowed.

#### Assigning A/G/U/C (function `assign_seq_for_pattern`)

For each S/W pattern, it enumerates all `2^N` assignments:

- At each W position: choose **A or U**
- At each S position: choose **G or C**

For each resulting sequence `s` (length 9), it applies filters:

1. **All four nucleotides must appear at least once**
- `s` must contain **A, G, U, and C** somewhere.

3. **Hairpin exclusion (for N ≥ 9, so yes for 9‑mers)**
- It looks at the would‑be stem formed by the outer 3 pairs:
- Pairs `s[0]–s[-1]`, `s[1]–s[-2]`, `s[2]–s[-3]`.

- Two scores:
- `score` using normal base pairs (AU, UA, GC, CG, GU, UG),
- `score_comp` using a slightly different “complement” rule (allows AC/CA).

- If **all three** outer positions look pairable in either sense (`score ≥ 3` or `score_comp ≥ 3`), the sequence is **discarded** – it’s too hairpin‑like.

5. **No internal “UUU” or mirrored “AAA”**
- For positions `i = 2 ... N-3` (internal):
- Reject if there is `U U U` at `s[i:i+3]`.
- Reject if there is `A A A` at the **mirrored** positions at the 5′ side.

- So we avoid long **internal U‑tracts and A‑tracts**, which are structurally problematic.

Sequences passing these tests are collected into a pool.

#### Final refinement filters (function `refine_seq`)

After pooling all candidates from all patterns, `refine_seq` applies four more rules:

1. `if6_kmry` **– avoid 6‑long degenerate runs**
- K = {G, U}
- M = {A, C}
- R = {A, G}
- Y = {C, U}
- It tracks the longest consecutive run for each group (K, M, R, Y).
   If **any group has length ≥ 6**, the sequence is **rejected**. → Roughly: no “6 bases in a row that all behave like purines, pyrimidines, etc.”

3. `if_ends_ww` **– both ends weak? reject**
- If both 5′ and 3′ ends are A/U (in any combination A/U vs A/U) the sequence is rejected.
- Ensures **at least one end has a G or C**, which stabilises the designed KL and reduces problematic 3′/5′ fluff.

5. `if_ends_3w` **– 3 weak bases at either end**
- If positions 0–2 or positions N-3–N-1 are all A/U, reject.

7. `if_ends_pair_u` **– more detailed end hairpin checks**
- For N ≥ 9, it looks for patterns like:
- 5′ A at the start or 3′ U at the end,
- short A‑rich segments at 5′ that could strongly pair with 3′,
- 3′ terminal U or UU that could form a little end‑stem with the 5′ side.

- If these patterns indicate likely terminal pairing, the sequence is rejected.



---

## Step 2 – `build_rna_conflict_graphV2.py`

**Goal:** from the raw pool, throw away self‑complementary and guide‑complementary sequences, then build a **graph of pairwise conflicts** based on strong complementarity.

### How to run it (9‑nt example)

Basic usage with defaults:


```
python build_rna_conflict_graphV2.py --num-nt 9
```
- `--num-nt 9` tells it we are working with 9‑mers.
- If you don’t give `--input`, it assumes: (the output of Step 1).
- Default guide(s): `DEFAULT_GUIDES = ["CACGAAGUCAAUAC"]`
   (this is the 15‑nt handle used in the bricks scripts). You can override:


```
python build_rna_conflict_graphV2.py \
--num-nt 9 \
--input RNAPool_9nt.txt \
--guide CACGAAGUCAAUAC --guide GGGAAAUUU \
```
The run time of the code could be up to 2 hours. 

With defaults you get:

- `FilteredPool_9nt_out.txt`
- `ConflictGraph_9nt_edges_out.txt`

### Step 2a – read and pre‑filter

1. **Read all sequences** (must have length `num_nt`).
2. **Remove palindromes / self‑complementary sequences**
- A sequence counts as “palindromic” if it is complementary to itself under the same rules used for pairwise complementarity (including GU wobble).

4. **Remove sequences that conflict with guide sequences**For each guide:
- Generate all contiguous 9‑mers from the guide.
- For each pool sequence, check if it is “too complementary” to **any** of those 9‑mers (with all reverse‑complement orientations).
- If yes → discard the pool sequence.

This ensures **no KL candidate binds strongly to the constant handle (or any other guide you specify).**

### Step 2b – complementarity rules and the conflict graph

Now you have a filtered list of sequences `seqs`. The script builds an undirected graph:

- **One node per sequence.**
- An **edge (i, j)** means sequences **i** and **j** are “too complementary” under the custom rules.

Under the hood:

- Base pairs allowed:
- Watson–Crick: **A‑U, U‑A, G‑C, C‑G**
- GU wobble: **G‑U, U‑G**

- For each pair `(s1, s2)` it checks **all four combinations**:
- `s1` vs `s2`
- `s1` vs reverse‑complement(s2)
- reverse‑complement(s1) vs `s2`
- reverse‑complement(s1) vs reverse‑complement(s2)

The **complementarity criteria** are exactly summarized in the header of the script:

Two sequences (or their RCs) are considered a **conflicting pair** if **any** of these holds:

1. **Global pairing**
- There exists an alignment with:
- at least **N–1 paired positions**, and
- **> 2 GC pairs**.

3. **Very long consecutive pairing**
- Alignment has a **longest consecutive run ≥ N–1**, and
- **> 1 GC pair** in that alignment.

5. **Medium consecutive pairing**
- Longest consecutive run ≥ **N–2**, and
- **> 2 GC pairs**.

7. **Shorter but GC‑rich**
- Longest consecutive run ≥ **N–3**, and
- **≥ 4 GC pairs**.

9. **Bulged pairing (N ≥ 9)**
- There exists a **single‑bulge alignment** of length **N–1** between the two sequences (one base skipped on one strand) with **≥ 3 GC pairs**.

For **N = 9** these thresholds are:

- N–1 = **8**, N–2 = **7**, N–3 = **6**.
- So two 9‑mers conflict if e.g.:
- They can form ~8 bp with several GC pairs, **or**
- They have a contiguous 8‑bp run with >1 GC, **or**
- They have a contiguous 7‑bp run with >2 GC, **or**
- They have a contiguous 6‑bp run with ≥4 GC, **or**
- They can form an **8‑bp single‑bulged** alignment with ≥3 GC.

This is more nuanced than just “no 5‑bp perfect complement”—it balances **length of interaction** and **GC richness**, which controls the expected stability.

**Output of Step 2 (for N=9):**

- **Filtered sequences** (post palindrome + guide filter)
   → `FilteredPool_9nt_out.txt`
- **Conflict graph edge list**
   → `ConflictGraph_9nt_edges_out.txt`(each line: `i<TAB>j`, with `0 ≤ i < j`)



---

## Step 3 – `select_from_conflict_graph.py`

**Goal:** choose a large subset of sequences that have **no edges between them** – a big independent set in the conflict graph.

### How to run it (9‑nt example)

If you used the defaults in Step 2:


```
python select_from_conflict_graph.py \
--seq-file   FilteredPool_9nt_out.txt \
--graph-file ConflictGraph_9nt_edges_out.txt \
--rounds 1000
```
It will write: `SelectedPool_from_graph_out.txt` 

Each line:


```
index<TAB>sequence
```
(the index matches the line number in `FilteredPool_9nt_out.txt`).

Key arguments:

- `--seq-file` : filtered sequences from Step 2
- `--graph-file` : edge list from Step 2
- `--rounds` : number of randomized greedy attempts (default 1000)
- `--seed` : random seed (optional but good for reproducibility)
- `--output` : base name for output file (default `SelectedPool_from_graph.txt` → `SelectedPool_from_graph_out.txt`)

### Step 3a – greedy independent set

Core routine: `greedy_independent_set(adj, rng, seed=None)`

Algorithm:

1. Start with:
- `remaining = {0, 1, ..., n-1}` (all nodes),
- `selected = []`.

3. Optionally add a **seed** set of already‑independent nodes.
4. While there are nodes remaining:
- Pick a random node `v` from `remaining`.
- Add `v` to `selected`.
- Remove `v` **and all its neighbors** from `remaining`.

This produces a **maximal** independent set (you can’t add any more nodes without creating conflicts), but not necessarily maximum.

In `main`, the script runs this many times:


```
for round_idx in range(args.rounds):
selected = greedy_independent_set(adj, rng)
if len(selected) > len(best_indices):
best_indices = selected
```
It keeps the largest independent set it has seen so far (`best_indices`).

### Step 3b – optimization pass (`optimize_select_graph`)

To push the size a bit further, there’s a second phase that tries to improve `best_indices`:

- Print “Before optimization: k” where `k = len(best_indices)`.
- If `k == 0`, nothing to do.
- Otherwise:
11365. Set `count_unchanged = 0`.
11366. Repeat while `count_unchanged ≤ 120 * current_size`:
- If current size ≤ 1, break.
- Create a **seed** by dropping the first node:
- Run `greedy_independent_set(adj, rng, seed=seed)`.
- If the new set is **larger**, accept it and reset `count_unchanged = 0`.
- Otherwise, accept it but increment `count_unchanged`.

- You can hit **Ctrl‑C** to stop earlier.

Intuition: “kick out” one sequence and re‑grow greedily: sometimes this lets you escape a local optimum and add more sequences overall.

**Final output:** a large set of **mutually non‑conflicting 9‑nt KL candidates** in `SelectedPool_from_graph_out.txt`. (Typically, >220 candidates; if you are lucky, >230.)



---

The result is a **large, algorithmically optimized set of 9‑nt KL sequences** that:

- behave nicely on their own,
- don’t bind strongly to the constant handle,
- and are predicted not to bind strongly to each other in any orientation.

### **Convert the selected pool by** `pool_rna_complement.py` **([Download link](https://www.dropbox.com/scl/fi/95uvzoltp6w05fdnmhy1y/pool_rna_complement.py?rlkey=n1vfwjkzu39k5zmxreq3a7ny3&dl=0)).**

This final pool is what you feed into tile‑design scripts (2D/3D bricks), which then assign one 9‑nt bulge (with its reverse‑complement loop) to each bKL edge in the structure.

The **tile‑design scripts** expect a **KL pool file** where each non‑empty line has **at least two** RNA sequences, and they treat the **first column** as the 9‑nt bulge sequence.

To get a compatible pool file from the graph‑selected list, we run:


```
python pool_rna_complement.py SelectedPool_from_graph_out.txt
```
- **Input format:** `index<TAB>sequence`
- **Output format:** `sequence<TAB>complement_sequence`

The output file is written next to the input, with `"_complement"` inserted before the extension.
