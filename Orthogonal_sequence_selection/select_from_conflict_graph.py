"""
select_from_conflict_graph.py


Script 2 of 2.


This script:
- Reads a list of node labels (e.g. RNA sequences), one per line.
- Reads an undirected conflict graph as an edge list "i<TAB>j" (0-based indices).
- Interprets this as a graph where an edge means "these two nodes conflict".
- Uses randomized greedy independent-set heuristics to find a large subset
  of nodes with no edges between them (a large independent set).
- Performs an additional optimization pass that repeatedly drops one node
  and recomputes a greedy independent set to try to enlarge the solution.
- Writes the selected nodes (index and label) to an output file.


This script is intentionally agnostic to the meaning of the nodes and edges,
so it can be reused for other applications.


Example
-------
    python select_from_conflict_graph.py \\
        --seq-file FilteredPool_9nt_out.txt \\
        --graph-file ConflictGraph_9nt_edges_out.txt \\
        --rounds 2000
"""


import argparse
import os
import random
import sys
from typing import List, Set, Tuple




def add_suffix_to_path(path: str, suffix: str) -> str:
    """
    Add a suffix before the file extension.


    "file.txt" + "_out" -> "file_out.txt"
    "file" + "_out" -> "file_out"
    """
    base, ext = os.path.splitext(path)
    if ext:
        return f"{base}{suffix}{ext}"
    return f"{path}{suffix}"




def read_labels(path: str) -> List[str]:
    """Read node labels (e.g. sequences), one per line."""
    labels: List[str] = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            labels.append(s)
    print(f"{len(labels)} node labels read from {path}.")
    return labels




def read_edge_list(path: str, num_nodes: int) -> List[Set[int]]:
    """
    Read an undirected edge list from a file.


    Each non-empty, non-comment line must be:
        i<TAB>j   or   i<SPACE>j
    where 0 <= i, j < num_nodes.


    Returns:
        adjacency list: list of sets of neighbors.
    """
    adj: List[Set[int]] = [set() for _ in range(num_nodes)]
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.replace(",", " ").split()
            if len(parts) != 2:
                raise ValueError(f"Invalid edge line {line!r} in {path!r}")
            i, j = int(parts[0]), int(parts[1])
            if i < 0 or j < 0 or i >= num_nodes or j >= num_nodes:
                raise ValueError(
                    f"Edge ({i}, {j}) out of range for num_nodes={num_nodes}"
                )
            if i == j:
                continue  # ignore self-loops if present
            adj[i].add(j)
            adj[j].add(i)
    print(f"Conflict graph read from {path}.")
    return adj




def greedy_independent_set(adj: List[Set[int]], rng: random.Random, seed: List[int] = None) -> List[int]:
    """
    Compute a greedy randomized independent set on the given conflict graph.


    Optionally starts from a given seed (which is assumed to already be
    independent). The algorithm then adds more vertices at random,
    always maintaining independence.


    Returns:
        A list of node indices forming an independent set.
    """
    n = len(adj)
    remaining = set(range(n))
    selected: List[int] = []


    if seed:
        for v in seed:
            if v in remaining:
                selected.append(v)
                remaining.remove(v)
                remaining.difference_update(adj[v]) #The difference_update() method removes the items that exist in both sets.


    while remaining:
        v = rng.choice(tuple(remaining))
        selected.append(v)
        remaining.remove(v)
        remaining.difference_update(adj[v])


    return selected




def optimize_select_graph(
    adj: List[Set[int]],
    current_selected: List[int],
    rng: random.Random,
) -> List[int]:
    """
    Graph-based analogue of the original OptimizeSelect.


    Repeatedly drops the first node from the current selection and
    recomputes a greedy randomized independent set using the remaining
    selection as a seed. Stops after ~120 * current_size non-improving
    iterations or on Ctrl+C.


    Returns:
        A (typically larger) list of node indices forming an independent set.
    """
    print("Optimizing selected set on conflict graph (Ctrl+C to stop early)...")
    print(f"Before optimization: {len(current_selected)}")


    if not current_selected:
        print("Selected set is empty; nothing to optimize.")
        return current_selected


    count_unchanged = 0
    try:
        while count_unchanged <= len(current_selected) * 120:
            if len(current_selected) <= 1:
                break
            seed = current_selected[1:]
            new_selected = greedy_independent_set(adj, rng, seed=seed)


            if len(new_selected) > len(current_selected):
                current_selected = new_selected
                print(f"* improved to size {len(current_selected)}")
                count_unchanged = 0
            else:
                current_selected = new_selected
                count_unchanged += 1
    except KeyboardInterrupt:
        print("\nOptimization interrupted by user.")


    print(f"After optimization: {len(current_selected)}")
    return current_selected




def write_selected(path: str, selected_indices: List[int], labels: List[str]) -> None:
    """
    Write selected nodes to file.


    Each line: "index<TAB>label".
    """
    with open(path, "w") as f:
        for idx in selected_indices:
            f.write(f"{idx}\t{labels[idx]}\n")




def parse_args(argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Select a large independent set given labels and a conflict graph."
    )
    parser.add_argument(
        "--seq-file",
        type=str,
        required=True,
        help="File with node labels (e.g. RNA sequences), one per line.",
    )
    parser.add_argument(
        "--graph-file",
        type=str,
        required=True,
        help="File with conflict edges, one 'i<tab>j' per line (0-based indices).",
    )
    parser.add_argument(
        "--rounds",
        type=int,
        default=1000,
        help="Number of randomized independent-set rounds to run before optimization. "
        "Default: 1000.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None, # Python seeds the generator from system entropy: Typically a mix of OS randomness (os.urandom) and current time.
        help="Random seed for reproducible selection. If not set, system randomness is used.",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Base name for the output file of selected nodes. "
        "Default: SelectedPool_from_graph.txt (the script will append '_out').",
    )
    return parser.parse_args(argv)




def main(argv: List[str]) -> None:
    args = parse_args(argv)
    labels = read_labels(args.seq_file)
    num_nodes = len(labels)


    if num_nodes == 0:
        print("No labels found; nothing to do.")
        return


    adj = read_edge_list(args.graph_file, num_nodes)


    if args.output is None:
        base_output = "SelectedPool_from_graph.txt"
    else:
        base_output = args.output
    output_path = add_suffix_to_path(base_output, "_out")


    rng = random.Random(args.seed)


    best_indices: List[int] = []


    print(f"Running {args.rounds} randomized graph rounds (Ctrl+C to stop early)...")
    try:
        for round_idx in range(args.rounds):
            if round_idx % 10 == 0:
                print(f"Round {round_idx + 1}")
            selected = greedy_independent_set(adj, rng)
            if len(selected) > len(best_indices):
                best_indices = selected
                print(
                    f"  New best set size: {len(best_indices)} "
                    f"(after round {round_idx + 1})"
                )
    except KeyboardInterrupt:
        print("\nGraph rounds interrupted by user.")


    if not best_indices:
        print("No independent set found; nothing to write.")
        return


    best_indices = optimize_select_graph(adj, best_indices, rng)
    write_selected(output_path, best_indices, labels)


    print(f"Final set size: {len(best_indices)}")
    print(f"Selected nodes written to: {output_path}")




if __name__ == "__main__":
    main(sys.argv[1:])
