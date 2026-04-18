import glob
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Parameters
L = 8
D = 3
N = L**D
TARGET_T = 0.2

# Same as always
BASE_DIR = os.path.dirname(__file__)
style_file = os.path.join(BASE_DIR, 'lib', 'science.mplstyle')
if os.path.exists(style_file):
    plt.style.use(style_file)

FIG_DIR = os.path.join('..', 'figures')
RES4_DIR = os.path.join('..', 'results', 'part4')
RES6_DIR = os.path.join('..', 'results', 'part6')
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RES6_DIR, exist_ok=True)

pt_files = sorted(glob.glob(os.path.join(RES4_DIR, 'timeseries_part4_s*.dat')))
met_files = sorted(glob.glob(os.path.join(RES6_DIR, 'timeseries_part6_s*.dat')))


def sample_id(path):
    m = re.search(r's(\d+)', os.path.basename(path))
    if m:
        return m.group(1)
    return os.path.basename(path)


def read_timeseries(path):
    return pd.read_csv(
        path,
        comment='#',
        sep=r'\s+',
        header=None,
        names=['mcs', 'k', 'T', 'Q', 'Eavg']
    )


def nearest_temperature(df, target):
    temps = np.sort(df['T'].unique())
    return float(temps[np.argmin(np.abs(temps - target))])


pt_map = {sample_id(f): f for f in pt_files}
met_map = {sample_id(f): f for f in met_files}

common_samples = sorted(set(pt_map) & set(met_map))
if not common_samples:
    raise RuntimeError('No common samples found between part4 and part6.')


pt_frames = []
met_frames = []
temp_rows = []

for sid in common_samples:
    pt_df = read_timeseries(pt_map[sid])
    used_T = nearest_temperature(pt_df, TARGET_T)
    pt_df = pt_df[np.isclose(pt_df['T'], used_T, atol=1e-8)].copy()
    pt_df['q2'] = (pt_df['Q'].astype(float) / float(N))**2
    pt_df['run_q2'] = np.cumsum(pt_df['q2']) / np.arange(1, len(pt_df) + 1)
    pt_df['sample'] = f's{sid}'
    pt_frames.append(pt_df[['sample', 'mcs', 'T', 'q2', 'run_q2']])

    met_df = read_timeseries(met_map[sid]).copy()
    met_df['q2'] = (met_df['Q'].astype(float) / float(N))**2
    met_df['run_q2'] = np.cumsum(met_df['q2']) / np.arange(1, len(met_df) + 1)
    met_df['sample'] = f's{sid}'
    met_frames.append(met_df[['sample', 'mcs', 'T', 'q2', 'run_q2']])

    temp_rows.append({
        'sample': f's{sid}',
        'target_T': TARGET_T,
        'used_T_pt': used_T,
        'used_T_met': float(met_df['T'].iloc[0])
    })


pt_all = pd.concat(pt_frames, ignore_index=True)
met_all = pd.concat(met_frames, ignore_index=True)

pt_avg = (
    pt_all.groupby('mcs', as_index=False)
    .agg(mean_run_q2=('run_q2', 'mean'))
    .sort_values('mcs')
)

met_avg = (
    met_all.groupby('mcs', as_index=False)
    .agg(mean_run_q2=('run_q2', 'mean'))
    .sort_values('mcs')
)

summary = pd.merge(
    pt_avg.rename(columns={'mean_run_q2': 'pt_mean_run_q2'}),
    met_avg.rename(columns={'mean_run_q2': 'met_mean_run_q2'}),
    on='mcs',
    how='outer'
).sort_values('mcs')

summary.to_csv(os.path.join(RES6_DIR, 'summary_part6_running_q2.dat'), sep=' ', index=False)
pd.DataFrame(temp_rows).to_csv(os.path.join(RES6_DIR, 'part6_temperature_map.dat'), sep=' ', index=False)


# Plotting
# Plot every number of samples to avoid overcrowding the plot
step_size = 1000  # Plot every 100th sample
plt.figure()

for sample, grp in pt_all.groupby('sample'):
    plt.plot(grp.iloc[::step_size]['mcs'], grp.iloc[::step_size]['run_q2'], color='0.75', linewidth=1.0, alpha=0.9)

for sample, grp in met_all.groupby('sample'):
    plt.plot(grp.iloc[::step_size]['mcs'], grp.iloc[::step_size]['run_q2'], color='tab:orange', linestyle='--', linewidth=1.0, alpha=0.45)

plt.plot([], [], color='0.75', linewidth=1.2, label='PT samples')
plt.plot([], [], color='tab:orange', linestyle='--', linewidth=1.2, label='Metropolis samples')

# Average also written every step_size:
plt.plot(pt_avg.iloc[::step_size]['mcs'], pt_avg.iloc[::step_size]['mean_run_q2'], color='black', linewidth=2.2, label='PT average')
plt.plot(met_avg.iloc[::step_size]['mcs'], met_avg.iloc[::step_size]['mean_run_q2'], color='tab:red', linestyle='--', linewidth=2.2, label='Metropolis average')

plt.xlabel('Monte Carlo steps')
plt.ylabel(r'Running average of $\langle q^2 \rangle$')
plt.title(r'Part 6: convergence of $\langle q^2 \rangle$ at low temperature')
plt.legend(loc='center right')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'part6_q2_convergence_T0p2.pdf'))
plt.close()


print('Written:')
print('  ../results/part6/summary_part6_running_q2.dat')
print('  ../results/part6/part6_temperature_map.dat')
print('  ../figures/part6_q2_convergence_T0p2.pdf')