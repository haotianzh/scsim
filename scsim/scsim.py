import os
import sys
import importlib.resources
import numpy as np
import subprocess as sp
from contextlib import contextmanager
from . import util
from .tree import BaseTree


@contextmanager
def get_executable_path(provided_path=None):
    """
    Context manager to resolve the binary path safely using importlib.
    Handles cross-platform naming (scsim.exe vs scsim).
    """
    if provided_path:
        yield provided_path
        return
    binary_name = "scsim.exe" if sys.platform == "win32" else "scsim"
    try:
        with importlib.resources.path("scsim.bin", binary_name) as bin_path:
            yield str(bin_path)
    except (ImportError, ModuleNotFoundError):
        local_path = os.path.join(os.path.dirname(__file__), "bin", binary_name)
        if os.path.exists(local_path):
            yield local_path
        else:
            raise FileNotFoundError(
                f"Could not locate {binary_name} in package resources or {local_path}"
            )


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
    tree = util.relabel(tree, offset=1)
    cn = []
    reads_wild = []
    reads_mut = []
    tg = []

    # Write the temporary tree file
    with open(tmpfile, "w") as out:
        out.write(tree.output(branch_length_func=lambda x: x.branch))

    # Resolve executable and run simulation
    with get_executable_path(executable) as exe_path:
        # Construct command as a list (safer than shell=True string)
        cmd = [
            str(exe_path),
            str(tmpfile),
            str(n_site),
            str(n_vaiant_per_site),
            str(error),
            str(dropout),
            str(doublet),
            str(rate_cn_gain),
            str(rate_cn_loss),
            str(recurrent),
            "0",  # Placeholder arg from original code
            str(dropout_cell_variance),
            str(coverage_mean),
            str(coverage_std),
            str(random_seed),
            str(int(beta_binomial)),
        ]

        # Run the subprocess
        res = sp.run(
            cmd,
            shell=False,  # set to False for cross-platform safety with list args
            stdout=sp.PIPE,
            stderr=sp.PIPE,  # Capture stderr just in case
            text=True,
            check=True,  # Raise exception if C++ binary fails
        )

    # Parse Output
    for line in res.stdout.splitlines():
        line = line.strip()
        if "Number of sites" in line:
            # Example line: "Number of sites: 100"
            try:
                nsite = int(line.split(":")[1])
            except IndexError:
                continue

        if "Read Count (0,1)" in line:
            # Example line: "Read Count (0,1) = (10, 5)"
            parts = line.split("=")[1].strip()[1:-1].split(",")
            reads_wild.append(int(parts[0]))
            reads_mut.append(int(parts[1]))

        if "Copy Number =" in line:
            # Example line: "Copy Number = 2"
            cn.append(int(line.split("=")[1].strip()))

        if "Genotype =" in line:
            # Example line: "Genotype = (1,0)"
            g = line.split("=")[1].strip()[1:-1]
            g0_str, g1_str = g.split(",")
            # Handle potential non-int genotypes if necessary, assuming int here
            tg.append([int(g0_str), int(g1_str)])

    # Reshape outputs
    # Ensure nsite was found to avoid reshape errors
    if "nsite" not in locals():
        raise RuntimeError("Simulation failed: 'Number of sites' not found in output.")

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
