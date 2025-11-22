import os
import numpy as np
import subprocess as sp
from . import util
from .tree import BaseTree


def get_executable():
    executable = os.path.join(os.path.dirname(__file__), "bin", "scsim")
    return executable


def load_tree(tree_file):
    with open(tree_file) as f:
        nwk = f.readline().strip()
    tree = util.from_newick(nwk)
    return tree


def simulate(
    tree,
    n_site,
    n_vaiant_per_site=1,
    error=0.01,
    dropout=0.2,
    dropout_cell_variance=0,
    coverage_mean=10,
    coverage_std=5,
    doublet=0,
    recurrent=0,
    rate_cn_gain=0.05,
    rate_cn_loss=0.01,
    beta_binomial=False,
    random_seed=42,
    tmpfile="tmp_tree.nwk",
    executable=None,
):
    if executable is None:
        executable = os.path.join(os.path.dirname(__file__), "bin", "scsim")
    assert os.path.exists(executable), "scism not found."
    tree = util.relabel(tree, offset=1)
    cn = []
    reads_wild = []
    reads_mut = []
    tg = []
    with open(tmpfile, "w") as out:
        out.write(tree.output(branch_length_func=lambda x: x.branch))
    res = sp.run(
        f"{executable} {tmpfile} {n_site} {n_vaiant_per_site} {error} {dropout} {doublet} {rate_cn_gain} {rate_cn_loss} {recurrent} 0 {dropout_cell_variance} {coverage_mean} {coverage_std} {random_seed} {int(beta_binomial)}",
        shell=True,
        stdout=sp.PIPE,
        text=True,
    )
    for line in res.stdout.splitlines():
        if "Number of sites" in line:
            nsite = int(line.strip().split(":")[1])
        if "Read Count (0,1)" in line:
            reads_wild.append(
                int(line.strip().split("=")[1].strip()[1:-1].split(",")[0])
            )
            reads_mut.append(
                int(line.strip().split("=")[1].strip()[1:-1].split(",")[1])
            )
        if "Copy Number =" in line:
            cn.append(int(line.strip().split("=")[1].strip()))
        if "Genotype =" in line:
            g = line.strip().split("=")[1].strip()[1:-1]
            g0, g1 = (g.split(",")[0]), int(g.split(",")[1])
            tg.append([g0, g1])
    reads_wild = np.array(reads_wild).reshape(nsite, -1)
    reads_mut = np.array(reads_mut).reshape(nsite, -1)
    cn = np.array(cn).reshape(nsite, -1)
    tg = np.array(tg, dtype=int).reshape(nsite, -1, 2)
    return {
        "wild_counts": reads_wild,
        "mutant_counts": reads_mut,
        "true_genotype": tg,
        "true_tree": tree.output(branch_length_func=lambda x: x.branch),
    }
