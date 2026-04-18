import glob
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Parameters
L = 8
D = 3
N = L**D
LAST_INTERVAL = 4
MCS_LO = 300000
MCS_HI = 1000001
TARGET_TEMPS = [0.2, 0.5, 1.0, 2.0]


# Same as always
BASE_DIR = os.path.dirname(__file__)
style_file = os.path.join(BASE_DIR, 'lib', 'science.mplstyle')
if os.path.exists(style_file):
    plt.style.use(style_file)

FIG_DIR = os.path.join('..', 'figures')
RES4_DIR = os.path.join('..', 'results', 'part4')
RES5_DIR = os.path.join('..', 'results', 'part5')
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RES5_DIR, exist_ok=True)


summary_files = sorted(glob.glob(os.path.join(RES4_DIR, 'summary_part4_s*.dat')))
ts_files = sorted(glob.glob(os.path.join(RES4_DIR, 'timeseries_part4_s*.dat')))


def sample_label(path):
    base = os.path.basename(path)
    return base.replace('summary_', '').replace('timeseries_', '').replace('.dat', '')


# Load summary data (last interval only)
summary_frames = []
for f in summary_files:
    df = pd.read_csv(
        f,
        comment='#',
        sep=r'\s+',
        header=None,
        names=['interval', 'k', 'T', 'ndata', 'mean_q2', 'err_q2', 'mean_e', 'err_e']
    )
    df = df[df['interval'] == LAST_INTERVAL].copy()
    df['sample'] = sample_label(f)
    summary_frames.append(df)

summary_df = pd.concat(summary_frames, ignore_index=True).sort_values(['sample', 'T'])

avg_df = (
    summary_df.groupby(['k', 'T'], as_index=False)
    .agg(
        n_samples=('sample', 'nunique'),
        mean_q2=('mean_q2', 'mean'),
        std_q2=('mean_q2', 'std'),
        mean_e=('mean_e', 'mean'),
        std_e=('mean_e', 'std'),
        mean_err_q2=('err_q2', 'mean'),
        mean_err_e=('err_e', 'mean')
    )
    .sort_values('T')
)

avg_df['sem_q2'] = avg_df['std_q2'].fillna(0.0) / np.sqrt(avg_df['n_samples'])
avg_df['sem_e']  = avg_df['std_e'].fillna(0.0)  / np.sqrt(avg_df['n_samples'])

avg_out = os.path.join(RES5_DIR, 'summary_part5_last_interval_avg.dat')
avg_df.to_csv(avg_out, sep=' ', index=False)


# Plot q^2 vs T 
plt.figure()
for sample, grp in summary_df.groupby('sample'):
    plt.plot(grp['T'], grp['mean_q2'], marker='o', alpha=0.45,
             linewidth=1.2, label=sample)
plt.errorbar(
    avg_df['T'], avg_df['mean_q2'], yerr=avg_df['sem_q2'],
    fmt='o-', color='black', linewidth=2.0, capsize=3, label='disorder average'
)
plt.xlabel('Temperature T')
plt.ylabel(r'$\langle q^2 \rangle$')
plt.title(r'Part 5: last-interval $\langle q^2 \rangle$ vs $T$')
plt.legend(loc='best')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'part5_q2_vs_T_last_interval.pdf'))
plt.close()


# Plot E/N vs T
plt.figure()
for sample, grp in summary_df.groupby('sample'):
    plt.plot(grp['T'], grp['mean_e'], marker='o', alpha=0.45,
             linewidth=1.2, label=sample)
plt.errorbar(
    avg_df['T'], avg_df['mean_e'], yerr=avg_df['sem_e'],
    fmt='o-', color='black', linewidth=2.0, capsize=3, label='disorder average'
)
plt.xlabel('Temperature T')
plt.ylabel(r'$\langle E \rangle / N$')
plt.title(r'Part 5: last-interval $\langle E \rangle / N$ vs $T$')
plt.legend(loc='best')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'part5_e_vs_T_last_interval.pdf'))
plt.close()


# Build temperature map 
ref_df = pd.read_csv(
    ts_files[0],
    comment='#',
    sep=r'\s+',
    header=None,
    names=['mcs', 'k', 'T', 'Q', 'Eavg']
)
ref_df = ref_df[(ref_df['mcs'] >= MCS_LO) & (ref_df['mcs'] < MCS_HI)].copy()
available_temps = np.sort(ref_df['T'].unique())

chosen_rows  = []
chosen_temps = []
for target in TARGET_TEMPS:
    actual = float(available_temps[np.argmin(np.abs(available_temps - target))])
    chosen_rows.append({'target_T': target, 'used_T': actual})
    chosen_temps.append(actual)

pd.DataFrame(chosen_rows).to_csv(
    os.path.join(RES5_DIR, 'part5_temperature_map.dat'),
    sep=' ', index=False
)


# P(q) histograms
colors = plt.cm.tab10(np.linspace(0, 1, max(10, len(ts_files))))

saved_hists = []
for target, actual in zip(TARGET_TEMPS, chosen_temps):
    fig, ax = plt.subplots()

    for i, ts_file in enumerate(ts_files):
        df = pd.read_csv(
            ts_file,
            comment='#',
            sep=r'\s+',
            header=None,
            names=['mcs', 'k', 'T', 'Q', 'Eavg']
        )
        df = df[(df['mcs'] >= MCS_LO) & (df['mcs'] < MCS_HI)]
        df = df[np.isclose(df['T'], actual, atol=1e-8)].copy()
        if df.empty:
            continue
        qnorm = df['Q'].to_numpy(dtype=float) / float(N)
        ax.hist(
            qnorm,
            bins=40,
            density=True,
            histtype='step',
            linewidth=1.6,
            color=colors[i],
            label=sample_label(ts_file)
        )

    ax.set_title(
        r'Part 5: $P(q)$, last interval $[3\times10^5,\,10^6]$ MCS'
        f'\n$T = {actual:.4f}$'
    )
    ax.set_xlabel(r'$q = Q/N$')
    ax.set_ylabel(r'$P(q)$')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.25)
    fig.tight_layout()

    fname = f'part5_Pq_T{target:g}.pdf'.replace('.', 'p').replace('pp', 'p.')
    fname = f'part5_Pq_T{str(target).replace(".", "p")}.pdf'
    out_path = os.path.join(FIG_DIR, fname)
    fig.savefig(out_path)
    plt.close(fig)
    saved_hists.append(out_path)


print('Written:')
print(f'  {avg_out}')
print('  ../results/part5/part5_temperature_map.dat')
print('  ../figures/part5_q2_vs_T_last_interval.pdf')
print('  ../figures/part5_e_vs_T_last_interval.pdf')
for p in saved_hists:
    print(f'  {p}')