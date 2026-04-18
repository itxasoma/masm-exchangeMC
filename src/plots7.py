import glob
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
style_file = os.path.join(BASE_DIR, 'lib/science.mplstyle')
if os.path.exists(style_file):
    plt.style.use(style_file)

RES = "../results/part7"
FIG = "../figures"
os.makedirs(FIG, exist_ok=True)


def load_autocorr(pattern):
    files = sorted(glob.glob(pattern))
    out = {}
    for f in files:
        tag = os.path.basename(f).replace('.dat', '')
        df = pd.read_csv(
            f, comment='#', sep=r'\s+', header=None,
            names=['lag', 'rho', 'tau_run']
        )
        out[tag] = df
    return out


def plot_rho_on_ax(ax, data_dict, color, label):
    lags = None
    all_rho = []
    for df in data_dict.values():
        ax.plot(df['lag'], df['rho'], color=color, alpha=0.20, linewidth=0.7)
        all_rho.append(df['rho'].values)
        lags = df['lag'].values
    if lags is not None and all_rho:
        ax.plot(lags, np.mean(all_rho, axis=0),
                color=color, linewidth=1.8, label=label)


def plot_taurun_on_ax(ax, data_dict, color, label):
    lags = None
    all_tau = []
    for df in data_dict.values():
        ax.plot(df['lag'], df['tau_run'], color=color, alpha=0.20, linewidth=0.7)
        all_tau.append(df['tau_run'].values)
        lags = df['lag'].values
    if lags is not None and all_tau:
        ax.plot(lags, np.mean(all_tau, axis=0),
                color=color, linewidth=1.8, label=label)


def final_tauint(data_dict):
    return np.array([df['tau_run'].iloc[-1] for df in data_dict.values()])


# Load
pt_q  = load_autocorr(f"{RES}/autocorr_q_pt_s*.dat")
pt_e  = load_autocorr(f"{RES}/autocorr_e_pt_s*.dat")
met_q = load_autocorr(f"{RES}/autocorr_q_met_s*.dat")
met_e = load_autocorr(f"{RES}/autocorr_e_met_s*.dat")

if not any([pt_q, pt_e, met_q, met_e]):
    raise FileNotFoundError(
        "No autocorrelation .dat files in ../results/part7/. "
        "Run  make autocorr7  first."
    )

COLOR_PT  = 'steelblue'
COLOR_MET = 'darkorange'

# Figure 1: rho(lag) for q 
fig, ax = plt.subplots()
if pt_q:
    plot_rho_on_ax(ax, pt_q,  COLOR_PT,  'PT (mean)')
if met_q:
    plot_rho_on_ax(ax, met_q, COLOR_MET, 'Metropolis (mean)')
ax.axhline(0, color='k', linewidth=0.6, linestyle='--')
ax.set_xlabel('Lag (MC steps)')
ax.set_ylabel(r'$\rho(t)$')
ax.set_title(r'Autocorrelation of $q$  ($T=0.2$)')
ax.legend()
ax.set_xlim(left=0)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{FIG}/part7_autocorr_q_T0p2.pdf")
plt.close()
print(f"Written: {FIG}/part7_autocorr_q_T0p2.pdf")

# Figure 2: rho(lag) for E 
fig, ax = plt.subplots()
if pt_e:
    plot_rho_on_ax(ax, pt_e,  COLOR_PT,  'PT (mean)')
if met_e:
    plot_rho_on_ax(ax, met_e, COLOR_MET, 'Metropolis (mean)')
ax.axhline(0, color='k', linewidth=0.6, linestyle='--')
ax.set_xlabel('Lag (MC steps)')
ax.set_ylabel(r'$\rho(t)$')
ax.set_title(r'Autocorrelation of $E/N$  ($T=0.2$)')
ax.legend()
ax.set_xlim(left=0)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{FIG}/part7_autocorr_E_T0p2.pdf")
plt.close()
print(f"Written: {FIG}/part7_autocorr_E_T0p2.pdf")

# Figure 3: running tau_int vs lag
fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))

for ax, pt_d, met_d, title in zip(
    axes,
    [pt_q,  pt_e],
    [met_q, met_e],
    [r'$q$', r'$E/N$']
):
    if pt_d:
        plot_taurun_on_ax(ax, pt_d,  COLOR_PT,  'PT')
    if met_d:
        plot_taurun_on_ax(ax, met_d, COLOR_MET, 'Metropolis')
    ax.set_xlabel('Lag (MC steps)')
    ax.set_ylabel(r'$\tau_{\mathrm{int}}$ (running sum)')
    ax.set_title(f'Observable: {title}')
    ax.legend()
    ax.grid(True, alpha=0.3)

plt.suptitle(r'Running $\tau_{\mathrm{int}}$ at $T=0.2$')
plt.tight_layout()
plt.savefig(f"{FIG}/part7_tauint_running.pdf")
plt.close()
print(f"Written: {FIG}/part7_tauint_running.pdf")

# Summary table 
rows = []
for method, dq, de in [('PT', pt_q, pt_e), ('Metropolis', met_q, met_e)]:
    tq = final_tauint(dq) if dq else np.array([float('nan')])
    te = final_tauint(de) if de else np.array([float('nan')])
    rows.append({
        'method':         method,
        'tau_int_q_mean': float(np.mean(tq)),
        'tau_int_q_std':  float(np.std(tq)),
        'tau_int_e_mean': float(np.mean(te)),
        'tau_int_e_std':  float(np.std(te)),
    })

summary = pd.DataFrame(rows)
out_csv = f"{RES}/tauint_summary.dat"
summary.to_csv(out_csv, sep=' ', index=False)
print(f"Written: {out_csv}")
print()
print(summary.to_string(index=False))