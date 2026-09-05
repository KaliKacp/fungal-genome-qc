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

if __name__ == "__main__":
    os.makedirs("results/plots", exist_ok=True)
    lengths = load_seq_lengths("data/raw/genome.fasta")
    plot_contig_length_distribution(lengths, "results/plots/contig_length_distribution.png")
    print(f"Wygenerowano wykres dla {len(lengths)} sekwencji.")