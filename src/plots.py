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
FERDI_FILE = os.path.join(RESULTS_DIR, 'ferdi.D')

L = 20
N = L * L
DISCARD_MCS = 1000


def load_mc_data(filename):
    data = np.loadtxt(filename, comments='#')
    mcs = data[:, 0].astype(int)
    k = data[:, 1].astype(int)
    T = data[:, 2]
    E = data[:, 4]
    return mcs, k, T, E


def mc_energy_vs_T(filename):
    mcs, k, T, E = load_mc_data(filename)

    keep = mcs >= DISCARD_MCS
    mcs = mcs[keep]
    k = k[keep]
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


def load_ferdi(filename):
    data = np.loadtxt(filename)
    T = data[:, 0]
    e = data[:, 1]
    return T, e


def plot_energy_comparison():
    Tmc, emc, d_emc = mc_energy_vs_T(TS_FILE)

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

    out_file = os.path.join(FIG_DIR, 'energy_vs_temperature_part2.pdf')
    plt.savefig(out_file)
    plt.close()
    print(f'Generated {os.path.basename(out_file)}')


if __name__ == '__main__':
    print('Generating plots...')
    plot_energy_comparison()
    print('Done!')
