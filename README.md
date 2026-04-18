# Exchange Monte Carlo (Parallel Tempering) — Ising Spin Glass

_(MASM — Mètodes Avançats de Simulació Molecular, 2025-2026)_
_(Itxaso Muñoz-Aldalur)_

Fortran implementation of the Parallel Tempering algorithm for the 3D Ising spin glass
H(S) = -sum_{<i,j>} J_ij S_i S_j, validated on the 2D ferromagnetic Ising model.

## Compile

    cd src
    make clean && make

This builds: exchange_mc (parts 1,2,4,5), exchange_mc_part3 (part 3),
exchange_mc_part6 (part 6), main_binning.x, main_binning_part4.x,
main_autocorr_part7.x.

## Run & plot (each part)

Part 1 (test):
```
   make run
```

Part 2 (2D Ising, and Ferdinand validation):

    gfortran -O2 -o ferdinand/ferdinand.x ferdinand/ferdinand_matteo.f 
    make run2 && make bin2 && make figs2

Parts 3-5 (3D spin glass, cluster):

    sbatch 3.run.sh   # 1000 disorder samples, results in results/part3/
    sbatch 4.run.sh   # 5 samples (SLURM array), results in results/part4/
    make bin4 && make figs4 && make figs5

Part 6 (Metropolis at T=0.2, cluster):

    sbatch 6.run.sh   # 5 samples (SLURM array), results in results/part6/
    make figs6

Part 7 (autocorrelation times):

    make part7        # compiles, runs Fortran autocorr on parts 4+6, plots

## Reproduce plots from stored results

Heavy timeseries are stored as .zip in results/part*/. Unzip before plotting:

    unzip results/part2/timeseries_part2.dat.zip -d results/part2/
    unzip results/part4/timeseries_part4_s0*.dat.zip -d results/part4/

Then: make figs2  /  make figs3  /  make figs4  /  make figs5  /  make figs6  /  make figs7

## Python environment

    python3 -m venv .venv && source .venv/bin/activate
    pip install -r src/lib/requirements.txt
