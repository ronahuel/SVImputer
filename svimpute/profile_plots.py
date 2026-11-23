"""Render descriptive figures from exact VCF counts."""
import numpy as np

from .profile import BINS


def draw_figures(output, classes, types, lengths):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator, PercentFormatter

    def canvas(title, width=7.6, height=4.6):
        fig, ax = plt.subplots(figsize=(width, height), constrained_layout=True)
        ax.set_title(title, loc="left", pad=14)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_axisbelow(True)
        ax.grid(axis="y", alpha=.2)
        return fig, ax

    def save(fig, name):
        fig.savefig(output / name, dpi=180)
        plt.close(fig)

    names = sorted(types, key=lambda name: (-types[name], name))
    fig, ax = canvas("Structural variants by type", width=max(7.6, .8 * len(names)))
    bars = ax.bar(names, [types[name] for name in names], color="#306a91")
    ax.bar_label(bars, padding=3)
    ax.margins(y=.18)
    ax.set_ylabel("SV ALT entries")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    if not names:
        ax.text(.5, .5, "No structural variants in this scope", ha="center", transform=ax.transAxes)
    save(fig, "sv_types.png")

    fig, ax = canvas("Structural variant lengths by type")
    for name in names:
        values = sorted(lengths[name])
        if not values:
            continue
        n = sum(lengths[name].values())
        fractions = np.cumsum([lengths[name][v] for v in values]) / n
        # Include the first vertical rise, including for singleton distributions.
        ax.step([values[0], *values], [0., *fractions], where="post",
                marker="." if len(values) <= 12 else None, label=f"{name} (n={n:,})")
    known = sum(sum(group.values()) for group in lengths.values())
    unknown = sum(types.values()) - known
    ax.set(xlabel=f"SV length (bp, log scale)\nKnown length: {known:,}    Unknown length: {unknown:,}",
           ylabel="Cumulative share within type", ylim=(-.025, 1.05))
    ax.yaxis.set_major_formatter(PercentFormatter(1.))
    ax.set_xscale("log")
    if known:
        ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1, .5))
    else:
        ax.set_xlim(1, 10)
        ax.text(.5, .5, "No SV lengths available", ha="center", transform=ax.transAxes)
    save(fig, "sv_lengths.png")

    fig, ax = canvas("Alleles per positional locus")
    centers = np.arange(len(BINS))
    width = .35
    for shift, scope, color in [(-.5, "All loci", "#306a91"), (.5, "Loci with SVs", "#d07744")]:
        bars = ax.bar(centers + shift*width, [classes[scope][label] for label in BINS],
                      width=width, label=scope, color=color)
        ax.bar_label(bars, padding=3)
    ax.set_xticks(centers, BINS)
    ax.set_ylabel("Loci")
    ax.set_ylim(0, max(1, max(classes["All loci"].values(), default=0) * 1.22))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
    save(fig, "alleles_per_locus.png")
