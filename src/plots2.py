# See environment requirements in README.md
# python3 -m venv .venv
# source .venv/bin/activate
# python3 -m pip install --upgrade pip
# python3 -m pip install -r lib/requirements.txt

# gfortran -O2 -o ferdinand.x ferdinand/ferdinand_matteo.f

import numpy as np
import matplotlib.pyplot as plt
import os

BASE_DIR = os.path.dirname(__file__)
RESULTS_DIR = os.path.join(BASE_DIR, '../results')
FIG_DIR = os.path.join(BASE_DIR, '../figures')

style_file = os.path.join(BASE_DIR, 'lib/science.mplstyle')
if os.path.exists(style_file):
    plt.style.use(style_file)

os.makedirs(FIG_DIR, exist_ok=True)

TS_FILE = os.path.join(RESULTS_DIR, 'timeseries_part2.dat')
FERDI_FILE  = os.path.join(BASE_DIR, '../results/ferdi.D')  # stays at root
SUMMARY_FILE = os.path.join(RESULTS_DIR, 'summary_part2.dat')
BINNING_FILE = os.path.join(RESULTS_DIR, 'binning_part2.dat')
SWAP_FILE = os.path.join(RESULTS_DIR, 'swap_stats_part2.dat')

L = 20
N = L * L
DISCARD_MCS = 10000
TS_PLOT_STRIDE = 100


def load_mc_data(filename):
    data = np.loadtxt(filename, comments='#')
    mcs = data[:, 0].astype(int)
    k = data[:, 1].astype(int)
    T = data[:, 2]
    Q = data[:, 3]
    E = data[:, 4]
    return mcs, k, T, Q, E


def mc_energy_vs_T_raw(filename):
    mcs, k, T, Q, E = load_mc_data(filename)

    keep = mcs >= DISCARD_MCS
    T = T[keep]
    E = E[keep]

    temps = np.unique(T)
    mean_e = []
    err_e = []

    for temp in temps:
        x = E[T == temp] / N
        mean_e.append(np.mean(x))
        if len(x) > 1:
            err_e.append(np.std(x, ddof=1) / np.sqrt(len(x)))
        else:
            err_e.append(0.0)

    return temps, np.array(mean_e), np.array(err_e)


def load_summary(filename):
    data = np.loadtxt(filename, comments='#')
    T = data[:, 1]
    mean_e = data[:, 4]
    err_e = data[:, 6]
    return T, mean_e, err_e


def load_summary_from_binning(filename, min_blocks=16, rel_tol=0.05):
    data = np.loadtxt(filename, comments='#')
    temps = np.unique(data[:, 1])

    T_out = []
    mean_out = []
    err_out = []

    for temp in temps:
        mask = np.abs(data[:, 1] - temp) < 1e-12
        rows = data[mask]

        rows = rows[np.argsort(rows[:, 2])]

        block_size = rows[:, 2]
        nblocks = rows[:, 3]
        mean_spin = rows[:, 5]
        err_spin = rows[:, 7]

        good = nblocks >= min_blocks
        if np.any(good):
            block_size = block_size[good]
            mean_spin = mean_spin[good]
            err_spin = err_spin[good]

        idx_plateau = len(err_spin) - 1

        if len(err_spin) >= 3:
            for i in range(1, len(err_spin) - 1):
                r1 = abs(err_spin[i] - err_spin[i - 1]) / max(abs(err_spin[i]), 1e-16)
                r2 = abs(err_spin[i + 1] - err_spin[i]) / max(abs(err_spin[i + 1]), 1e-16)
                if r1 < rel_tol and r2 < rel_tol:
                    idx_plateau = i
                    break

        T_out.append(temp)
        mean_out.append(mean_spin[idx_plateau])
        err_out.append(err_spin[idx_plateau])

    return np.array(T_out), np.array(mean_out), np.array(err_out)


def load_ferdi(filename):
    data = np.loadtxt(filename)
    T = data[:, 0]
    e = data[:, 1]
    return T, e


def choose_plot_temps(temps, targets=(2.0, 2.2, 2.6)):
    chosen = []
    for t in targets:
        tsel = temps[np.argmin(np.abs(temps - t))]
        if not any(abs(tsel - x) < 1e-12 for x in chosen):
            chosen.append(tsel)
    return np.array(chosen)


def plot_energy_comparison():
    if os.path.exists(BINNING_FILE):
        Tmc, emc, d_emc = load_summary_from_binning(BINNING_FILE)
    elif os.path.exists(SUMMARY_FILE):
        Tmc, emc, d_emc = load_summary(SUMMARY_FILE)
    else:
        Tmc, emc, d_emc = mc_energy_vs_T_raw(TS_FILE)

    plt.figure()
    plt.errorbar(Tmc, emc, yerr=d_emc, fmt='o', capsize=3, label='Exchange MC')

    if os.path.exists(FERDI_FILE):
        Tf, ef = load_ferdi(FERDI_FILE)
        mask = (Tf >= Tmc.min() - 1e-12) & (Tf <= Tmc.max() + 1e-12)
        plt.plot(Tf[mask], ef[mask], '-', lw=1.8, label='Ferdinand-Fisher')

    plt.xlabel('Temperature')
    plt.ylabel('Energy per spin')
    plt.title('2D Ising: MC vs exact')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    out_file = os.path.join(FIG_DIR, 'part2_energy_vs_temperature.pdf')
    plt.savefig(out_file)
    plt.close()
    print(f'Generated {os.path.basename(out_file)}')


def plot_binning_curves():
    if not os.path.exists(BINNING_FILE):
        return

    data = np.loadtxt(BINNING_FILE, comments='#')
    temps = np.unique(data[:, 1])
    chosen = choose_plot_temps(temps)

    plt.figure()
    for temp in chosen:
        mask = np.abs(data[:, 1] - temp) < 1e-12
        block_size = data[mask, 2]
        err_spin = data[mask, 7]
        plt.plot(block_size, err_spin, 'o-', label=f'T={temp:.2f}')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Block size')
    plt.ylabel('Binned error of energy per spin')
    plt.title('Binning analysis')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    out_file = os.path.join(FIG_DIR, 'part2_binning_energy.pdf')
    plt.savefig(out_file)
    plt.close()
    print(f'Generated {os.path.basename(out_file)}')


def plot_energy_timeseries():
    if not os.path.exists(TS_FILE):
        return

    mcs, k, T, Q, E = load_mc_data(TS_FILE)
    keep = mcs >= DISCARD_MCS
    mcs = mcs[keep]
    T = T[keep]
    E = E[keep]

    temps = np.unique(T)
    chosen = choose_plot_temps(temps)

    plt.figure()
    for temp in chosen:
        mask = np.abs(T - temp) < 1e-12
        plt.plot(
            mcs[mask][::TS_PLOT_STRIDE],
            (E[mask] / N)[::TS_PLOT_STRIDE],
            label=f'T={temp:.2f}',
            alpha=0.9,
            lw=0.6
        )

    plt.xlabel('MCS')
    plt.ylabel('Energy per spin')
    plt.title('Energy time series')
    plt.grid(True, alpha=0.3)
    plt.legend(frameon=True)
    plt.tight_layout()

    out_file = os.path.join(FIG_DIR, 'part2_energy_timeseries.pdf')
    plt.savefig(out_file)
    plt.close()
    print(f'Generated {os.path.basename(out_file)}')


def plot_overlap_timeseries():
    if not os.path.exists(TS_FILE):
        return

    mcs, k, T, Q, E = load_mc_data(TS_FILE)
    keep = mcs >= DISCARD_MCS
    mcs = mcs[keep]
    T = T[keep]
    Q = Q[keep]

    temps = np.unique(T)
    chosen = choose_plot_temps(temps)

    plt.figure()
    for temp in chosen:
        mask = np.abs(T - temp) < 1e-12
        plt.plot(
            mcs[mask][::TS_PLOT_STRIDE],
            (Q[mask] / N)[::TS_PLOT_STRIDE],
            label=f'T={temp:.2f}',
            alpha=0.9,
            lw=0.6
        )

    plt.xlabel('MCS')
    plt.ylabel('q = Q/N')
    plt.title('Overlap time series')
    plt.grid(True, alpha=0.3)
    plt.legend(frameon=True)
    plt.tight_layout()

    out_file = os.path.join(FIG_DIR, 'part2_overlap_timeseries.pdf')
    plt.savefig(out_file)
    plt.close()
    print(f'Generated {os.path.basename(out_file)}')


def plot_swap_rates():
    if not os.path.exists(SWAP_FILE):
        return

    data = np.loadtxt(SWAP_FILE, comments='#')
    Tmid = 0.5 * (data[:, 1] + data[:, 2])
    rateA = data[:, 5]
    rateB = data[:, 8]

    plt.figure()
    plt.plot(Tmid, rateA, 'o-', label='Family A')
    plt.plot(Tmid, rateB, 's-', label='Family B')
    plt.axhline(0.3, color='k', ls='--', lw=1.0, alpha=0.7)
    plt.xlabel('Pair midpoint temperature')
    plt.ylabel('Swap acceptance rate')
    plt.title('Parallel tempering swap rates')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    out_file = os.path.join(FIG_DIR, 'part2_swap_rates.pdf')
    plt.savefig(out_file)
    plt.close()
    print(f'Generated {os.path.basename(out_file)}')


def print_energy_table():
    if os.path.exists(BINNING_FILE):
        Tmc, emc, d_emc = load_summary_from_binning(BINNING_FILE)
    elif os.path.exists(SUMMARY_FILE):
        Tmc, emc, d_emc = load_summary(SUMMARY_FILE)
    else:
        Tmc, emc, d_emc = mc_energy_vs_T_raw(TS_FILE)

    print('\nEnergy per spin table')
    print('T        E_MC/N        err(E/N)      E_exact        |diff|')
    print('-' * 62)

    if os.path.exists(FERDI_FILE):
        Tf, ef = load_ferdi(FERDI_FILE)
        for T, e_mc, de_mc in zip(Tmc, emc, d_emc):
            idx = np.argmin(np.abs(Tf - T))
            e_ex = ef[idx]
            diff = abs(e_mc - e_ex)
            print(f'{T:5.2f}   {e_mc: .10f}   {de_mc: .10f}   {e_ex: .10f}   {diff: .10f}')
    else:
        for T, e_mc, de_mc in zip(Tmc, emc, d_emc):
            print(f'{T:5.2f}   {e_mc: .10f}   {de_mc: .10f}')


if __name__ == '__main__':
    print('Generating plots...')
    plot_energy_comparison()
    plot_binning_curves()
    plot_energy_timeseries()
    plot_overlap_timeseries()
    plot_swap_rates()
    print_energy_table()
    print('Done!')