# scsim: Simulator of SNVs and CNAs in Single-Cell Data

**scsim** is a C++ simulator for generating single-cell sequencing data with both **single nucleotide variants (SNVs)** and **copy number alterations (CNAs)**, including both **copy gains** and **copy losses**.  
It is designed for testing and benchmarking phylogenetic and clonal inference methods under realistic single-cell mutation processes.

scsim was originally developed during the development of ScisTree (https://github.com/yufengwudcs/ScisTree). If you use scsim, please cite:


>*Accurate and Efficient Cell Lineage Tree Inference from Noisy Single Cell Data: the Maximum Likelihood Perfect Phylogeny Approach, Yufeng Wu, Bioinformatics, https://doi.org/10.1093/bioinformatics/btz676, 36(3):742-750, 2020.*

---

## 🔍 Overview

The simulator models the joint evolution of SNVs and CNAs along a given lineage tree.  
Starting from a diploid wild-type ancestor, mutations accumulate along each branch according to user-specified rates.

### Simulation process

1. **Initialization:**  
   The root genotype is initialized as `(2, 0)` (diploid homozygous wild-type).

2. **Evolution along the tree:**  
   - **Point Mutations (PMs):** Occur according to a Poisson process, with a user-specified mutation rate.  
     Each site can mutate only once, following the Infinite Sites (IS) assumption.  
   - **Copy Number Mutations (CNMs):** Randomly select an allele and either duplicate or delete it with equal probability.  
     Multiple CNMs can occur independently along the same branch.

3. **Read simulation:**  
   - For each extant (leaf) cell, allelic dropout is simulated by discarding alleles with a fixed dropout rate.  
   - If an allele is not dropped out, reads are generated with a small chance of sequencing error.  
   - Coverage and read-depth variance are adjustable parameters.

---

## ⚙️ Installation

Compile using any standard C++ compiler such as `g++`:

```bash
g++ src/scsim.cpp -O3 -o scsim
```
A compiled binary `scsim` is also available for Linux users. 

## 🚀 Usage

Run the simulator with the following command:

```bash
./scsim tree-file num-sites num-variants-per-site error-rate dropout-rate doublet \
cnv-rate-inc cnv-rate-dec rec-rate binary droput-cell-var ave-read-depth \
std-read-depth rnd-seed fbetabinomial
```

| Parameter            | Description                                                                                  |
| -------------------- | -------------------------------------------------------------------------------------------- |
| `tree-file`          | Input tree file in **Newick** format with branch lengths                                     |
| `num-sites`           | Number of somatic single nucleotide variants (SNVs)                                          |
| `num-variants-per-site` | Number of variants per site                                                                  |
| `error-rate`         | Genotype error rate                                                                          |
| `dropout-rate`       | Genotype dropout rate                                                                        |
| `doublet`            | Number of doublet errors (the total number of output genotypes will decrease by this amount) |
| `cnv-rate-inc`       | Rate of copy number increase events                                                          |
| `cnv-rate-dec`       | Rate of copy number decrease events                                                          |
| `rec-rate`           | Rate of recurrent mutations                                                                  |
| `binary`             | Output mode: `0` for binary, `1` for ternary (optional)                                      |
| `droput-cell-var`    | Dropout cell variation (optional)                                                            |
| `ave-read-depth`     | Average read depth (optional)                                                                |
| `std-read-depth`     | Standard deviation of read depth (optional)                                                  |
| `rnd-seed`           | Random seed (integer; optional)                                                              |
| `fbetabinomial`      | Use beta-binomial for read counts if set to `1`; `0` uses normal distribution (optional)     |



## 🧩 Example

Below is an example command demonstrating a typical simulation run:

```bash
./scsim tree-example.nwk 5 1 0.01 0.4 0 0.5 0.2 0 0 0 10 4 101013 0
```
Example Description:
This example simulates a dataset with: 10 cells, 5 somatic SNV sites (1 variant per site). Genotyping error rate of 1% and dropout rate of 40%. Copy number increase and decrease rates of 0.5 and 0.2, respectively. Mean sequencing depth of 10 reads per locus with 4 as standard deviation. Deterministic seed (101013) for reproducibility.
The output will contain simulated single-cell genotypes and read counts consistent with the provided tree and mutation parameters.



## Python Interface

We also provide a Python wrapper for the simulator. You can install the `scsim` package directly from PyPI.  
Note that a working **C++ compiler (`g++`) is required**. If your system does not have one, we recommend installing it via Conda:

```bash
conda install -c conda-forge cxx-compiler
```

Then, install `scsim` in your conda env:
```bash
pip install scsim
```

### Usage in Python
```python
import scsim 
tree = scsim.get_random_binary_tree(10)   # generate a tree with 10 leaves.
data = scsim.simulate(tree, n_site=100)   # simulate reads for 100 sites with default parameters.
```
