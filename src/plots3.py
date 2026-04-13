# So now, to plot the histogram, read the histogram files:
import glob, os
import numpy as np
import matplotlib.pyplot as plt

RESULTS_DIR = '../results/part3'
FIG_DIR     = '../figures'
L, D = 4, 3
N = L**D
NBINS = 2*N + 1
TARGET_TEMPS = [0.2, 0.5]

FILES = sorted(glob.glob(os.path.join(RESULTS_DIR, 'histogram_part3_s*.dat')))
print(f'Found {len(FILES)} histogram files')

# Accumulate counts over all disorder samples
# Data columns: k  T  Q  count  nmeas_total
all_data = []
for f in FILES:
    all_data.append(np.loadtxt(f, comments='#'))

data = np.vstack(all_data)   # shape: (nfiles * NT * NBINS, 5)

for T_target in TARGET_TEMPS:
    # Find the closest temperature
    temps = np.unique(data[:,1])
    T_use = temps[np.argmin(np.abs(temps - T_target))]

    mask  = np.abs(data[:,1] - T_use) < 1e-9
    sub   = data[mask]            # rows for this temperature

    Q_vals = sub[:,2].astype(int)
    counts = sub[:,3].astype(int)

    # Sum counts for each Q value across all samples
    Q_range = np.arange(-N, N+1)
    total_counts = np.zeros(len(Q_range), dtype=float)
    for q, c in zip(Q_vals, counts):
        total_counts[q + N] += c

    # Normalize to get P(q): convert Q -> q=Q/N, normalize so integral=1
    q_centers = Q_range / float(N)
    dq = 1.0 / N
    total = total_counts.sum() * dq
    Pq = total_counts / total

    plt.figure()
    plt.step(q_centers, Pq, where='mid', lw=1.8)
    plt.xlabel(r'$q = Q/N$')
    plt.ylabel(r'$P(q)$')
    plt.yscale('log')
    plt.title(f'3D spin glass, L={L}, T={T_use:.2f} ({len(FILES)} samples)')
    plt.xlim(0, 1)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, f'pq_T{T_target:.1f}.pdf'))
    plt.close()
    print(f'T={T_use:.3f}: saved pq_T{T_target:.1f}.pdf')



# Old script:

# import glob
# import os
# import numpy as np
# import matplotlib.pyplot as plt

# BASE_DIR = os.path.dirname(__file__)
# RESULTS_DIR = os.path.join(BASE_DIR, '../results/part3')
# FIG_DIR = os.path.join(BASE_DIR, '../figures')

# style_file = os.path.join(BASE_DIR, 'lib/science.mplstyle')
# if os.path.exists(style_file):
#     plt.style.use(style_file)

# os.makedirs(FIG_DIR, exist_ok=True)

# L = 4
# D = 3
# N = L ** D
# DISCARD_MCS = 3000
# TARGET_TEMPS = [0.2, 0.5]
# NBINS = 2 * N + 1

# FILES = sorted(glob.glob(os.path.join(RESULTS_DIR, 'timeseries_part3_s*.dat')))


# def load_one_file(filename):
#     data = np.loadtxt(filename, comments='#')
#     mcs = data[:, 0].astype(int)
#     T = data[:, 2]
#     Q = data[:, 3]
#     return mcs, T, Q


# def collect_q_values(files, target_temp, discard_mcs=3000):
#     all_q = []
#     used_files = 0
#     chosen_temp = None

#     for fname in files:
#         mcs, T, Q = load_one_file(fname)

#         temps_in_file = np.unique(T)
#         tsel = temps_in_file[np.argmin(np.abs(temps_in_file - target_temp))]
#         if chosen_temp is None:
#             chosen_temp = tsel

#         mask = (np.abs(T - tsel) < 1e-12) & (mcs >= discard_mcs)
#         if np.any(mask):
#             q = Q[mask] / float(N)
#             all_q.append(q)
#             used_files += 1

#     if not all_q:
#         return None, None, 0

#     return np.concatenate(all_q), chosen_temp, used_files


# def histogram_pq(q_values, nbins=NBINS):
#     edges = np.linspace(-1.0 - 1.0 / N, 1.0 + 1.0 / N, nbins + 1)
#     hist, edges = np.histogram(q_values, bins=edges, density=True)
#     centers = 0.5 * (edges[:-1] + edges[1:])
#     return centers, hist


# def save_histogram_data(temp_target, temp_used, q_values):
#     centers, hist = histogram_pq(q_values)
#     outname = os.path.join(FIG_DIR, f'pq_hist_T{temp_target:.1f}.dat')
#     arr = np.column_stack([centers, hist])
#     header = f'q_center Pq   target_T={temp_target:.3f} used_T={temp_used:.6f}'
#     np.savetxt(outname, arr, header=header)
#     print(f'Saved {os.path.basename(outname)}')
#     return centers, hist


# def plot_single_hist(temp_target, temp_used, q_values, nfiles):
#     centers, hist = save_histogram_data(temp_target, temp_used, q_values)

#     plt.figure()
#     plt.step(centers, hist, where='mid', lw=1.8)
#     plt.xlabel(r'$q = Q/N$')
#     plt.xlim(0, 1.0)
#     plt.yscale('log')
#     plt.ylabel(r'$P(q)$')
#     plt.title(f'3D spin glass, L=4, T={temp_used:.2f}')
#     plt.grid(True, alpha=0.3)
#     plt.tight_layout()

#     outname = os.path.join(FIG_DIR, f'pq_hist_T{temp_target:.1f}.pdf')
#     plt.savefig(outname)
#     plt.close()
#     print(f'Generated {os.path.basename(outname)} using {nfiles} samples')


# def plot_both_histograms(files):
#     plt.figure()

#     for temp_target in TARGET_TEMPS:
#         q_values, temp_used, nfiles = collect_q_values(files, temp_target, DISCARD_MCS)
#         if q_values is None:
#             print(f'No data found for target T={temp_target}')
#             continue

#         centers, hist = histogram_pq(q_values)
#         plt.step(centers, hist, where='mid', lw=1.8, label=f'T={temp_used:.2f} ({nfiles} samples)')

#     plt.xlabel(r'$q = Q/N$')
#     plt.ylabel(r'$P(q)$')
#     plt.xlim(0, 1.0)
#     plt.yscale('log')
#     plt.title('3D spin glass overlap histograms')
#     plt.grid(True, alpha=0.3)
#     plt.legend()
#     plt.tight_layout()

#     outname = os.path.join(FIG_DIR, 'pq_hist_compare_T02_T05.pdf')
#     plt.savefig(outname)
#     plt.close()
#     print(f'Generated {os.path.basename(outname)}')


# def main():
#     if not FILES:
#         raise FileNotFoundError(f'No files found in {RESULTS_DIR}')

#     print(f'Found {len(FILES)} time-series files')

#     for temp_target in TARGET_TEMPS:
#         q_values, temp_used, nfiles = collect_q_values(FILES, temp_target, DISCARD_MCS)
#         if q_values is None:
#             print(f'No data found for target T={temp_target}')
#             continue
#         plot_single_hist(temp_target, temp_used, q_values, nfiles)

#     plot_both_histograms(FILES)


# if __name__ == '__main__':
#     main()
