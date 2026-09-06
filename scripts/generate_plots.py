"""
Generuje wykresy diagnostyczne assembly na podstawie surowego pliku FASTA.
Wyjście: pliki PNG do folderu results/plots/, użyte później w raporcie HTML.
"""
import matplotlib
matplotlib.use("Agg")  # backend bez GUI - wymagane w środowisku headless (Codespace/CI)
import matplotlib.pyplot as plt
from Bio import SeqIO
import os

def load_seq_lengths(fasta_path: str) -> list[int]:
    return [len(record.seq) for record in SeqIO.parse(fasta_path, "fasta")]

def plot_contig_length_distribution(lengths: list[int], outpath: str):
    fig, ax = plt.subplots(figsize=(8, 5))
    sorted_lengths = sorted(lengths, reverse=True)
    ax.bar(range(len(sorted_lengths)), sorted_lengths, color="#2c7fb8")
    ax.set_yscale("log")
    ax.set_xlabel("Sekwencja (posortowana malejąco)")
    ax.set_ylabel("Długość (bp, skala log)")
    ax.set_title("Rozkład długości sekwencji w assembly H. erinaceus")
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)

def gc_content(seq) -> float:
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    return 100 * gc / len(seq) if len(seq) > 0 else 0.0

def plot_gc_per_sequence(records, outpath: str):
    names = [r.id for r in records]
    gc_values = [gc_content(str(r.seq)) for r in records]
    lengths = [len(r.seq) for r in records]

    genome_avg_gc = sum(l * g for l, g in zip(lengths, gc_values)) / sum(lengths)

    fig, ax = plt.subplots(figsize=(9, 5))
    colors = ["#d95f0e" if abs(g - genome_avg_gc) > 5 else "#2c7fb8" for g in gc_values]
    ax.bar(range(len(names)), gc_values, color=colors)
    ax.axhline(genome_avg_gc, color="black", linestyle="--", linewidth=1,
               label=f"Średnia genomu: {genome_avg_gc:.1f}%")
    ax.set_xlabel("Sekwencja")
    ax.set_ylabel("GC content (%)")
    ax.set_title("GC content per sekwencja — kontrola kontaminacji")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)

def plot_busco_summary(busco: dict, outpath: str):
    categories = ["Single-copy", "Duplicated", "Fragmented", "Missing"]
    values = [
        busco["single_copy_pct"],
        busco["duplicated_pct"],
        busco["fragmented_pct"],
        busco["missing_pct"],
    ]
    colors = ["#2c7fb8", "#7fcdbb", "#f4a261", "#d95f02"]

    fig, ax = plt.subplots(figsize=(9, 2.5))
    left = 0
    for cat, val, color in zip(categories, values, colors):
        ax.barh(0, val, left=left, color=color, label=f"{cat} ({val}%)")
        left += val

    ax.set_xlim(0, 100)
    ax.set_yticks([])
    ax.set_xlabel("Procent BUSCO groups (%)")
    ax.set_title(f"BUSCO Completeness — {busco['lineage']} (n={busco['total_buscos']})")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.4), ncol=4, frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)

if __name__ == "__main__":
    import json

    os.makedirs("results/plots", exist_ok=True)

    records = list(SeqIO.parse("data/raw/genome.fasta", "fasta"))
    lengths = [len(r.seq) for r in records]

    plot_contig_length_distribution(lengths, "results/plots/contig_length_distribution.png")
    plot_gc_per_sequence(records, "results/plots/gc_content_per_sequence.png")

    with open("results/summary.json") as f:
        summary = json.load(f)
    plot_busco_summary(summary["busco"], "results/plots/busco_summary.png")

    print(f"Wygenerowano wykresy dla {len(records)} sekwencji + wykres BUSCO.")