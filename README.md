# Contaminationination Detection with Wheeler Graphs

## Overview

This repository introduces a reference-guided approach for contamination detection in genomic datasets. Contrasting with traditional reference-based methods that typically rely on searching for specific genes across multiple, separate data structures derived from reference FASTA files, our method employs the creation of a Wheeler graph to enable efficient and accurate contamination detection.

---

## Prepare the Environment

Before running any of the code, make sure to

From the root directory, run `source wheeler_graph_env/bin/activate`

This command should have you set up for running any of the scripts in this repository except for those in `./Baselines` (see the README in there for instructions).

---

## Running

For running individual scripts, `cd` into the desired subfolder and take a look at the README.md file there for specific instructions.

For testing the entire pipeline, FASTA files -> Wheeler Graph -> Querying FASTQ file -> Results, run `test.py` from the root directory (this can take up to a couple of minutes).

---

## Credit

The build de bruijn graph code is built upon this amazing wheeler graph repo: [https://github.com/Kuanhao-Chao/Wheeler_Graph_Toolkit](https://github.com/Kuanhao-Chao/Wheeler_Graph_Toolkit)

Check it out to learn more about recognizing and generating all sorts of wheeler graphs!
