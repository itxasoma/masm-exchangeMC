# Exchange Monte Carlo (Parallel Tempering) for the Ising Spin Glass

_(for the subject: Advanced Simulation Methods / Monte Carlo Methods)_ 

_(by: itxasoma (Itxaso Muñoz-Aldalur))_ 

Fortran implementation of the Exchange Monte Carlo (Parallel Tempering) algorithm for the Ising spin glass with Hamiltonian $H(S) = -\sum_{\langle i,j\rangle} J_{ij} S_i S_j$, using periodic boundary conditions and Metropolis updates inside each replica. The project first requires validation on the 2D ferromagnetic Ising model by setting $J_{ij}=1$, and then production simulations for the 3D Gaussian spin glass using overlap observables and disorder averaging. 

The simulation runs two independent replica families, usually labeled $a$ and $b$, each with $N_T$ temperatures between $T_{\min}$ and $T_{\max}$, and attempts exchanges between neighboring temperatures every `nsw` Monte Carlo sweeps. This is the standard parallel tempering idea: replicas move through temperature space so low-temperature configurations can escape metastable valleys more efficiently than in fixed-temperature Metropolis sampling. 

## How to
1. **Go to the** ```src``` **directory:**
   ```bash
   cd src
   ```
2. **Compile & run the simulation with plots:**
   ```bash
   make clean
   make
   ```
    2.1. (Optional) initial test to check that the code is working:
    ```bash
    make run
    ```
    2.2. 2D Ising model simulation: run, binning and plots:
   
    To run localy, tune the nMCS in inputs/part2.in. Also run `ferdinand_matteo.f`for comparison
   
   ```bash
   make run2
   make binning
   gfortran -O2 -o ferdinand.x ferdinand/ferdinand_matteo.f 
   make figures2
   ```
   
    To run in the cluster (used pirineus3.csuc.cat):
   ```bash
   sbatch 2.run.sh
    ```
    
    2.3. 3D Ising model simulation
   
    To run in the cluster (used pirineus3.csuc.cat):
   ```bash
   sbatch 3.run.sh
    ```
    inputs/part3_samples/ and logs/part3 will be created containing the simulation's inputs, as well as results/part3 with the swap statistics and timeseries files.

3. **Install Python in a virtual environment:**

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install numpy
python3 -m pip install matplotlib
python3 -m pip install pandas
```

## Structure

```text
.
├── figures
│   ├── binning_energy_part2.pdf
│   ├── energy_timeseries_part2.pdf
│   ├── energy_vs_temperature_part2.pdf
│   ├── overlap_timeseries_part2.pdf
│   ├── pq_hist_compare_T02_T05_part3.pdf
│   ├── pq_hist_T0.2_part3.pdf
│   ├── pq_hist_T0.2.dat
│   ├── pq_hist_T0.5_part3.pdf
│   ├── pq_hist_T0.5.dat
│   └── swap_rates_part2.pdf
├── inputs
│   ├── part1.in
│   ├── part2.in
│   ├── part3.in
│   └── temps.dat
├── LICENSE
├── README.md
├── results
│   ├── binning_part2.dat
│   ├── ferdi.D
│   ├── summary_part2.dat
│   ├── swap_stats_part1.dat
│   ├── swap_stats_part2.dat
│   ├── swap_stats.dat
│   ├── timeseries_exchMC.dat
│   ├── timeseries_part1.dat
│   └── timeseries_part2.dat
├── src
│   ├── 2.run.sh
│   ├── 3.run.sh
│   ├── binning.f90
│   ├── bonds.f90
│   ├── exchangeMC.f90
│   ├── ferdinand
│   │   ├── constants.par
│   │   ├── ferdinand_matteo.f
│   │   └── implicit.sta
│   ├── lattice.f90
│   ├── lib
│   │   ├── requirements.txt
│   │   └── science.mplstyle
│   ├── main_binning.f90
│   ├── main.f90
│   ├── Makefile
│   ├── parameters.f90
│   ├── plots2.py
│   ├── plots3.py
│   ├── r1279
│   │   ├── gauss.f
│   │   ├── r1279.f90
│   │   ├── r1279block.h
│   │   └── ran2.f
│   └── rng_wrapper.f90
└──
