import pandas as pd
import matplotlib.pyplot as plt

quick = pd.read_csv("quick_standard_results.csv")
intro = pd.read_csv("intro_results.csv")
typeNames = {
    0: "Random",
    1: "ReverseSorted",
    2: "AlmostSorted"
}
for t in [0, 1, 2]:
    sub = quick[quick["type"] == t]
    plt.figure(figsize=(10, 6))
    plt.plot(sub["n"], sub["mean_us"])
    plt.title(f"Quick Sort (mean time) — {typeNames[t]}")
    plt.xlabel("Размер массива n")
    plt.ylabel("Время, микросекунды")
    plt.grid(True)
    plt.savefig(f"{typeNames[t]}Quick.png")

for t in [0, 1, 2]:
    sub = intro[intro["type"] == t]
    plt.figure(figsize=(10, 6))
    plt.plot(sub["n"], sub["mean_us"])
    plt.title(f"Intro Sort (mean time) — {typeNames[t]}")
    plt.xlabel("Размер массива n")
    plt.ylabel("Время, микросекунды")
    plt.grid(True)
    plt.savefig(f"{typeNames[t]}Intro.png")
