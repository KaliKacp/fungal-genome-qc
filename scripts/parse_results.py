"""
Parsuje wyniki seqkit i QUAST do jednej wspólnej struktury danych (dict),
którą później wykorzysta generator raportu HTML.
"""
import pandas as pd
import json
import sys

def parse_seqkit(path: str) -> dict:
    df = pd.read_csv(path, sep="\t")
    row = df.iloc[0]  # jeden plik wejściowy = jeden wiersz
    return {
        "num_seqs": int(row["num_seqs"]),
        "total_length": int(row["sum_len"]),
        "gc_percent": float(row["GC(%)"]),
        "n50": int(row["N50"]),
        "min_len": int(row["min_len"]),
        "max_len": int(row["max_len"]),
    }

def parse_quast(path: str) -> dict:
    df = pd.read_csv(path, sep="\t", index_col=0)
    col = df.columns[0]  # nazwa kolumny = nazwa assembly ("genome")
    return {
        "n_contigs": int(df.loc["# contigs", col]),
        "largest_contig": int(df.loc["Largest contig", col]),
        "n50": int(df.loc["N50", col]),
        "l50": int(df.loc["L50", col]),
        "n_per_100kbp": float(df.loc["# N's per 100 kbp", col]),
    }

def parse_busco(path: str) -> dict:
    with open(path) as f:
        data = json.load(f)
    results = data["results"]
    return {
        "lineage": data["lineage_dataset"]["name"],
        "complete_pct": float(results["Complete percentage"]),
        "single_copy_pct": float(results["Single copy percentage"]),
        "duplicated_pct": float(results["Multi copy percentage"]),
        "fragmented_pct": float(results["Fragmented percentage"]),
        "missing_pct": float(results["Missing percentage"]),
        "total_buscos": int(results["n_markers"]),
    }

if __name__ == "__main__":
    summary = {
        "seqkit": parse_seqkit("results/seqkit/stats.tsv"),
        "quast": parse_quast("results/quast/report.tsv"),
        "busco": parse_busco("results/busco/busco_hericium/run_fungi_odb10/short_summary.json"),
    }
    with open("results/summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))