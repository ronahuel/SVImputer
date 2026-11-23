"""Aggregate signed correlations, giving each locus equal weight."""
from collections import defaultdict
import csv

import numpy as np


def writer(stream, fields):
    table = csv.writer(stream, delimiter="\t", lineterminator="\n")
    table.writerow(fields)
    return table


class Accuracy:
    def __init__(self):
        self.cells = defaultdict(lambda: dict(total=0, n_available=0, n_paired=0,
                                               sum_available=0., sum_paired=0.,
                                               min_used=None, max_used=0))

    def add_locus(self, chrom, correlations, sizes):
        # Keys within each window are (repeat, ALT). Pair within a locus first.
        common = set.intersection(*(set(values) for values in correlations.values()))
        for k, values in correlations.items():
            for scope in (("genome", ""), ("chromosome", chrom)):
                cell = self.cells[scope, k]
                cell["total"] += 1
                used = sizes[k]
                cell["min_used"] = used if cell["min_used"] is None else min(cell["min_used"], used)
                cell["max_used"] = max(cell["max_used"], used)
                if values:
                    cell["n_available"] += 1
                    cell["sum_available"] += float(np.mean(list(values.values())))
                if common:
                    cell["n_paired"] += 1
                    cell["sum_paired"] += float(np.mean([values[key] for key in sorted(common)]))

    def rows(self):
        rows = []
        for ((scope, chrom), k), cell in sorted(self.cells.items()):
            n_paired, n_available = cell["n_paired"], cell["n_available"]
            rows.append(dict(scope=scope, chrom=chrom, k=k, min_used=cell["min_used"],
                             max_used=cell["max_used"], loci=cell["total"],
                             evaluable_loci=n_available, paired_loci=n_paired,
                             mean_r_available=cell["sum_available"] / n_available if n_available else None,
                             mean_r_paired=cell["sum_paired"] / n_paired if n_paired else None, best=False))
        groups = defaultdict(list)
        for row in rows:
            if row["mean_r_paired"] is not None:
                groups[(row["scope"], row["chrom"])].append(row)
        for group in groups.values():
            peak = max(r["mean_r_paired"] for r in group)
            for row in group:
                row["best"] = abs(row["mean_r_paired"] - peak) <= 1e-12
        return rows


def plot(rows, destination):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    genome = [r for r in rows if r["scope"] == "genome"]
    fig, ax = plt.subplots(figsize=(7.4, 4.6), constrained_layout=True)
    indexes = list(range(len(genome)))
    ax.plot(indexes, [np.nan if r["mean_r_paired"] is None else r["mean_r_paired"] for r in genome],
            marker="o", color="#306a91", linewidth=2)
    best = [i for i, r in enumerate(genome) if r["best"]]
    if best:
        ax.scatter(best, [genome[i]["mean_r_paired"] for i in best], s=78,
                   color="#b44b20", zorder=3, label="Highest mean")
        ax.legend(frameon=False)
        ax.set_title(f"Mean dosage accuracy across {genome[0]['paired_loci']} comparable loci", loc="left")
    else:
        ax.text(.5, .55, "No comparable loci with defined Pearson r",
                ha="center", transform=ax.transAxes)
    ax.set_xticks(indexes)
    ax.set_xticklabels([str(r["k"]) for r in genome], rotation=45 if len(genome) > 10 else 0)
    ax.set(xlabel="Maximum number of nearby LD SNPs", ylabel="Mean signed Pearson r", ylim=(-1.04, 1.04))
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.grid(axis="y", alpha=.22)
    fig.savefig(destination, dpi=180)
    plt.close(fig)
