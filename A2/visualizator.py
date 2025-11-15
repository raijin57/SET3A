import pandas as pd
import matplotlib.pyplot as plt

merge = pd.read_csv("merge_standard_results.csv")
hybrid = pd.read_csv("merge_insertion_results5.csv")
typeNames = {
    0: "Random",
    1: "ReverseSorted",
    2: "AlmostSorted"
}
for t in [0, 1, 2]:
    sub = merge[merge["type"] == t]
    plt.figure(figsize=(10, 6))
    plt.plot(sub["n"], sub["mean_us"])
    plt.title(f"Merge Sort (mean time) — {typeNames[t]}")
    plt.xlabel("Размер массива n")
    plt.ylabel("Время, микросекунды")
    plt.grid(True)
    plt.savefig(f"{typeNames[t]}Merge.png")

for t in [0, 1, 2]:
    sub = hybrid[hybrid["type"] == t]
    plt.figure(figsize=(10, 6))
    plt.plot(sub["n"], sub["mean_us"])
    plt.title(f"Merge+Insertion Sort (mean time) — {typeNames[t]}")
    plt.xlabel("Размер массива n")
    plt.ylabel("Время, микросекунды")
    plt.grid(True)
    plt.savefig(f"{typeNames[t]}Hybrid5.png")

hybrid = pd.read_csv("merge_insertion_results10.csv")
for t in [0, 1, 2]:
    sub = hybrid[hybrid["type"] == t]
    plt.figure(figsize=(10, 6))
    plt.plot(sub["n"], sub["mean_us"])
    plt.title(f"Merge+Insertion Sort (mean time) — {typeNames[t]}")
    plt.xlabel("Размер массива n")
    plt.ylabel("Время, микросекунды")
    plt.grid(True)
    plt.savefig(f"{typeNames[t]}Hybrid10.png")

hybrid = pd.read_csv("merge_insertion_results20.csv")
for t in [0, 1, 2]:
    sub = hybrid[hybrid["type"] == t]
    plt.figure(figsize=(10, 6))
    plt.plot(sub["n"], sub["mean_us"])
    plt.title(f"Merge+Insertion Sort (mean time) — {typeNames[t]}")
    plt.xlabel("Размер массива n")
    plt.ylabel("Время, микросекунды")
    plt.grid(True)
    plt.savefig(f"{typeNames[t]}Hybrid20.png")

hybrid = pd.read_csv("merge_insertion_results30.csv")
for t in [0, 1, 2]:
    sub = hybrid[hybrid["type"] == t]
    plt.figure(figsize=(10, 6))
    plt.plot(sub["n"], sub["mean_us"])
    plt.title(f"Merge+Insertion Sort (mean time) — {typeNames[t]}")
    plt.xlabel("Размер массива n")
    plt.ylabel("Время, микросекунды")
    plt.grid(True)
    plt.savefig(f"{typeNames[t]}Hybrid30.png")

hybrid = pd.read_csv("merge_insertion_results50.csv")
for t in [0, 1, 2]:
    sub = hybrid[hybrid["type"] == t]
    plt.figure(figsize=(10, 6))
    plt.plot(sub["n"], sub["mean_us"])
    plt.title(f"Merge+Insertion Sort (mean time) — {typeNames[t]}")
    plt.xlabel("Размер массива n")
    plt.ylabel("Время, микросекунды")
    plt.grid(True)
    plt.savefig(f"{typeNames[t]}Hybrid50.png")
