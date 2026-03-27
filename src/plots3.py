import glob
import os
import numpy as np
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(__file__)
RESULTS_DIR = os.path.join(BASE_DIR, '../results/part3')
FIG_DIR = os.path.join(BASE_DIR, '../figures')

style_file = os.path.join(BASE_DIR, 'lib/science.mplstyle')
if os.path.exists(style_file):
    plt.style.use(style_file)

os.makedirs(FIG_DIR, exist_ok=True)

L = 4
D = 3
N = L ** D
DISCARD_MCS = 3000
TARGET_TEMPS = [0.2, 0.5]
NBINS = 2 * N + 1

FILES = sorted(glob.glob(os.path.join(RESULTS_DIR, 'timeseries_part3_s*.dat')))


def load_one_file(filename):
    data = np.loadtxt(filename, comments='#')
    mcs = data[:, 0].astype(int)
    T = data[:, 2]
    Q = data[:, 3]
    return mcs, T, Q


def collect_q_values(files, target_temp, discard_mcs=3000):
    all_q = []
    used_files = 0
    chosen_temp = None

    for fname in files:
        mcs, T, Q = load_one_file(fname)

        temps_in_file = np.unique(T)
        tsel = temps_in_file[np.argmin(np.abs(temps_in_file - target_temp))]
        if chosen_temp is None:
            chosen_temp = tsel

        mask = (np.abs(T - tsel) < 1e-12) & (mcs >= discard_mcs)
        if np.any(mask):
            q = Q[mask] / float(N)
            all_q.append(q)
            used_files += 1

    if not all_q:
        return None, None, 0

    return np.concatenate(all_q), chosen_temp, used_files


def histogram_pq(q_values, nbins=NBINS):
    edges = np.linspace(-1.0 - 1.0 / N, 1.0 + 1.0 / N, nbins + 1)
    hist, edges = np.histogram(q_values, bins=edges, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, hist


def save_histogram_data(temp_target, temp_used, q_values):
    centers, hist = histogram_pq(q_values)
    outname = os.path.join(FIG_DIR, f'pq_hist_T{temp_target:.1f}.dat')
    arr = np.column_stack([centers, hist])
    header = f'q_center Pq   target_T={temp_target:.3f} used_T={temp_used:.6f}'
    np.savetxt(outname, arr, header=header)
    print(f'Saved {os.path.basename(outname)}')
    return centers, hist


def plot_single_hist(temp_target, temp_used, q_values, nfiles):
    centers, hist = save_histogram_data(temp_target, temp_used, q_values)

    plt.figure()
    plt.step(centers, hist, where='mid', lw=1.8)
    plt.xlabel(r'$q = Q/N$')
    plt.xlim(0, 1.0)
    plt.yscale('log')
    plt.ylabel(r'$P(q)$')
    plt.title(f'3D spin glass, L=4, T={temp_used:.2f}')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    outname = os.path.join(FIG_DIR, f'pq_hist_T{temp_target:.1f}.pdf')
    plt.savefig(outname)
    plt.close()
    print(f'Generated {os.path.basename(outname)} using {nfiles} samples')


def plot_both_histograms(files):
    plt.figure()

    for temp_target in TARGET_TEMPS:
        q_values, temp_used, nfiles = collect_q_values(files, temp_target, DISCARD_MCS)
        if q_values is None:
            print(f'No data found for target T={temp_target}')
            continue

        centers, hist = histogram_pq(q_values)
        plt.step(centers, hist, where='mid', lw=1.8, label=f'T={temp_used:.2f} ({nfiles} samples)')

    plt.xlabel(r'$q = Q/N$')
    plt.ylabel(r'$P(q)$')
    plt.xlim(0, 1.0)
    plt.yscale('log')
    plt.title('3D spin glass overlap histograms')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    outname = os.path.join(FIG_DIR, 'pq_hist_compare_T02_T05.pdf')
    plt.savefig(outname)
    plt.close()
    print(f'Generated {os.path.basename(outname)}')


def main():
    if not FILES:
        raise FileNotFoundError(f'No files found in {RESULTS_DIR}')

    print(f'Found {len(FILES)} time-series files')

    for temp_target in TARGET_TEMPS:
        q_values, temp_used, nfiles = collect_q_values(FILES, temp_target, DISCARD_MCS)
        if q_values is None:
            print(f'No data found for target T={temp_target}')
            continue
        plot_single_hist(temp_target, temp_used, q_values, nfiles)

    plot_both_histograms(FILES)


if __name__ == '__main__':
    main()
