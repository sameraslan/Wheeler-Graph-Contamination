# Querying a Wheeler Graph for Contamination

---

## Overview

This folder contains a script for querying a Wheeler graph represented in a `.dot` file format. The primary functionality of this tool is to search for specific nucleotide patterns within the graph and identify the corresponding organism IDs. Example `.dot` files are located in the `graphs` subfolder. Currently, the functions in read_and_query.py are used in our test function so there's no need to run anything here, but you can call read_and_query.py if you like and it'll look for the example pattern in the ./graphs/example.dot wheeler graph.

---

## Prerequisites

If you are in the wheeler_graph_env then you don't need to worry about this. Otherwise, you need:

- NetworkX library, installable via `pip install networkx`.
- Pydot library (for reading `.dot` files), installable via `pip install pydot`.

---

## Features

- **Pattern Querying**: Search for specific nucleotide patterns within the Wheeler graph.
- **Organism ID Extraction**: Identify and list organism IDs that correspond to the queried pattern.
- **FASTQ File Validation**: Ensures the integrity and format of the input FASTQ files.

---

## Example Usage

```python
# To query a specific pattern in a graph run this in the if __name__ == "__main__":
G = nx.nx_pydot.read_dot(DOTFILE)
pattern = "CTAGG"
organism_ids = query_graph(G, pattern)
print(f"Organism(s) with this Pattern: {organism_ids}")
```
