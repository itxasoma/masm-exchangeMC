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
    ```sbatch 2.run.sh
    ```
    2.3. 3D Ising model simulation
    To run in the cluster (used pirineus3.csuc.cat):
    ```sbatch 3.run.sh
    ```


## Structure

```text
.
├── src/
│   ├── rng_module.f90           # Random number generator interface
│   ├── parameters.f90           # Global parameters and inputs
│   ├── lattice.f90              # Lattice geometry and periodic boundary conditions
│   ├── bonds.f90                # Quenched couplings Jij (Gaussian or Jij = 1 test mode)
│   ├── energy.f90               # Energy, local field, and ΔE routines
│   ├── metropolis.f90           # Single-temperature Metropolis sweeps
│   ├── tempering.f90            # Replica exchange / temperature swap moves
│   ├── observables.f90          # Q, q, q^2, energy, acceptance statistics
│   ├── binning.f90              # Binning analysis for statistical errors
│   ├── io.f90                   # Output files and logging
│   ├── main.f90                 # Main program: Parallel Tempering simulation
│   ├── plots.ipynb              # Analysis and plots
│   └── Makefile
├── inputs/
│   ├── test_2d_ising.in         # Validation case: d = 2, Jij = 1
│   ├── sg_3d_L4.in              # 3D spin glass, small-size production
│   └── sg_3d_L8.in              # 3D spin glass, long equilibration run
├── out/                         # Output files created at runtime
├── plots/                       # Figures created from analysis scripts/notebooks
└── README.md
