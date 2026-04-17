import glob
import os
import pandas as pd
import matplotlib.pyplot as plt

summary_files = sorted(glob.glob("../results/part4/summary_part4_s*.dat"))

if not summary_files:
    raise FileNotFoundError("No summary_part4_s*.dat files found in ../results/part4/")

BASE_DIR = os.path.dirname(__file__)
style_file = os.path.join(BASE_DIR, 'lib/science.mplstyle')
if os.path.exists(style_file):
    plt.style.use(style_file)


frames = []
for f in summary_files:
    df = pd.read_csv(
    f,
    comment="#",
    sep=r"\s+",
    header=None,
    names=["interval", "k", "T", "ndata", "mean_q2", "err_q2", "mean_e", "err_e"]
    )   
    df["sample"] = os.path.basename(f)
    frames.append(df)

all_df = pd.concat(frames, ignore_index=True)

avg = (
    all_df.groupby(["interval", "k", "T"], as_index=False)
    .agg(
        mean_q2=("mean_q2", "mean"),
        err_q2=("err_q2", "mean"),
        mean_e=("mean_e", "mean"),
        err_e=("err_e", "mean"),
    )
    .sort_values(["interval", "T"])
)

interval_labels = {
    1: r"$10^4 \leq t < 3\times10^4$",
    2: r"$3\times10^4 \leq t < 10^5$",
    3: r"$10^5 \leq t < 3\times10^5$",
    4: r"$3\times10^5 \leq t < 10^6$",
}

plt.figure(figsize=(8,6))
for interval, grp in avg.groupby("interval"):
    plt.plot(grp["T"], grp["err_q2"], marker="o", label=interval_labels.get(interval, f"interval {interval}"))
plt.xlabel("Temperature T")
plt.ylabel(r"Error of $q^2$")
plt.title(r"Part 4: disorder-averaged binning error of $q^2$")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("../figures/part4_err_q2_vs_T.pdf")

plt.figure(figsize=(8,6))
for interval, grp in avg.groupby("interval"):
    plt.plot(grp["T"], grp["err_e"], marker="o", label=interval_labels.get(interval, f"interval {interval}"))
plt.xlabel("Temperature T")
plt.ylabel(r"Error of $E/N$")
plt.title(r"Part 4: disorder-averaged binning error of $E/N$")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("../figures/part4_err_e_vs_T.pdf")

avg.to_csv("../results/part4/summary_part4_avg.dat", sep=" ", index=False)
print("Written:")
print("  ../results/part4/summary_part4_avg.dat")
print("  ../figures/part4_err_q2_vs_T.pdf")
print("  ../figures/part4_err_e_vs_T.pdf")
