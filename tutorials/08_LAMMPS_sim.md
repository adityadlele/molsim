
# Tutorial 8: Getting Started with LAMMPS and Molecular Dynamics

Welcome to the world of Molecular Dynamics (MD)! In this tutorial, we will transition from writing our own Python loops to using professional, high-performance MD software. We will learn the basics of setting up a simulation, defining the physics, and running it on a supercomputer.

---

## 1. What is LAMMPS?

**LAMMPS** stands for **L**arge-scale **A**tomic/**M**olecular **M**assively **P**arallel **S**imulator. 

* It is an open-source MD code originally developed by Sandia National Laboratories.
* It is one of the most widely used simulation engines in the world, particularly tailored for **materials science**, soft matter, polymers, and solid-state physics. 
* *Note:* While LAMMPS *can* perform biological simulations, bio-related research (like simulating proteins, DNA, or lipid bilayers) typically relies on specialized codes such as **GROMACS** or **NAMD**, which are highly optimized for those specific molecular structures.

---

## 2. What Do You Need to Run a LAMMPS Simulation?

To successfully run an MD job, you need four core components:

1. **A Geometry File (or Generation Commands):** You need to tell the software where the atoms are located in 3D space. You can either read in a pre-made coordinate file (e.g., a `.data` file built in a molecule editor) or use LAMMPS commands to generate a crystal lattice or a random gas box directly in the input script.
2. **An MD Potential (Force Field):** You must define how the atoms interact with one another. This dictates the underlying physics of the simulation.
3. **Simulation Parameters:** You must specify the thermodynamic ensemble (NVE, NVT, NPT), the timestep, the target temperature/pressure, and the total duration of the simulation.
4. **A Job Submission Script:** Because MD requires heavy computation, you will not run this on your laptop. You will submit it to a supercomputer cluster (like Purdue's Anvil) using a scheduler script.



---

## 3. Our First System: Argon Gas

For our first simulation, we are going to model Argon gas. 

* **Historical Significance:** The very first realistic molecular dynamics simulation was performed by Aneesur Rahman in 1964 using liquid Argon. It is considered the "Hello World" of statistical mechanics and computational chemistry.
* **The Potential:** Argon is a noble gas, meaning it does not form chemical bonds. Its atomic interactions are perfectly described by the purely non-bonded **Lennard-Jones (LJ) potential**, which models both the van der Waals attraction at long distances and the Pauli repulsion at short distances.



* **The Ensemble:** We will run an **NVE** simulation (Microcanonical Ensemble). This means the **N**umber of particles, the **V**olume of the box, and the Total **E**nergy are strictly conserved. We will not use a thermostat; the atoms will simply move according to Newton's equations of motion.

---

## 4. The LAMMPS Input File Structure

LAMMPS scripts are traditionally named with an `in.` prefix (e.g., `in.argon`). The script is read top-to-bottom and is generally divided into four distinct sections: Initialization, System Definition, Settings, and Run.

Create a file named `in.argon` on your computer and copy the following code into it:

```lammps
# -------------------------------------------------------------------
# 1. INITIALIZATION
# -------------------------------------------------------------------
units           real            # time=fs, energy=kcal/mol, distance=Angstroms
dimension       3               # 3D simulation
boundary        p p p           # Periodic boundary conditions in x, y, z
atom_style      atomic          # Simplest style: atoms have mass, no charge/bonds

# -------------------------------------------------------------------
# 2. SYSTEM DEFINITION
# -------------------------------------------------------------------
# Create a cubic box of 50x50x50 Angstroms
region          mybox block 0 50 0 50 0 50
create_box      1 mybox         # Create a simulation box with 1 atom type

# Populate the box with 500 atoms randomly
create_atoms    1 random 500 12345 mybox 

# Define the mass of our atom type 1 (Argon = 39.95 g/mol)
mass            1 39.95

# -------------------------------------------------------------------
# 3. SIMULATION SETTINGS (The Physics)
# -------------------------------------------------------------------
# Define the Lennard-Jones potential and a cutoff of 10.0 Angstroms
pair_style      lj/cut 10.0

# pair_coeff: atom_type1 atom_type2 epsilon sigma
pair_coeff      1 1 0.238 3.405  

# Give atoms an initial velocity corresponding to 300 K
# 87287 is a random seed for the velocity generator
velocity        all create 300.0 87287 dist gaussian

# -------------------------------------------------------------------
# 4. RUN DEFINITION
# -------------------------------------------------------------------
# Apply the NVE integrator to all atoms
fix             1 all nve

# Tell LAMMPS what to print to the screen/log file every 100 steps
# step = timestep, temp = temperature, pe = potential energy, ke = kinetic energy, etotal = total energy
thermo          100
thermo_style    custom step temp pe ke etotal press

# Dump the atomic coordinates to a file every 500 steps so we can visualize it later
dump            1 all custom 500 dump.argon id type x y z

# Define the timestep (1.0 femtoseconds is standard for real units)
timestep        1.0

# Run the simulation for 10,000 steps
run             10000

```

---

## 5. Job Submission Script for Anvil

Supercomputers like Anvil use a workload manager called **Slurm** to queue and distribute jobs fairly among thousands of users. You must write a shell script that requests computing resources (like CPU cores and time) and tells the computer how to launch LAMMPS.

Create a file named `submit_lammps.sh` and copy the following code. **Make sure to update the account allocation code to match our class!**

```bash
#!/bin/bash
#SBATCH --job-name=argon_nve          # Name of your job
#SBATCH --account=your_allocation     # FIXME: Replace with your specific class allocation code
#SBATCH --partition=shared            # Use the shared queue
#SBATCH --nodes=1                     # Number of nodes requested
#SBATCH --ntasks=4                    # Number of CPU cores requested (4 is plenty for 500 atoms)
#SBATCH --time=00:10:00               # Max walltime (10 minutes)
#SBATCH --output=lammps_job_%j.out    # Standard output file (%j will be replaced by job ID)
#SBATCH --error=lammps_job_%j.err     # Standard error file

# 1. Load the required LAMMPS module on Anvil
module purge
module load lammps

# 2. Move to the directory where you submitted the job
cd $SLURM_SUBMIT_DIR

# 3. Run LAMMPS using the MPI runner (srun)
srun lmp -in in.argon

```

### How to Run Your Simulation

1. **Log in** to the Anvil cluster via SSH or the Open OnDemand web portal.
2. **Upload** both your `in.argon` and `submit_lammps.sh` files to a new folder in your home or scratch directory.
3. Open a terminal in that folder and **submit the job** by typing:
`sbatch submit_lammps.sh`
4. You can **check the status** of your job in the queue using:
`squeue -u $USER`
5. **Review the output:** Once finished, LAMMPS will generate a `log.lammps` file containing your thermodynamic data (Temperature, Energy, Pressure) and a `dump.argon` file containing the 3D trajectory of your atoms for visualization.


## 6. Moving to Realistic Ensembles: NVT and NPT

While the NVE ensemble is fantastic for verifying that your physics are correct (since total energy should be perfectly conserved), real-world laboratory experiments are rarely isolated from their surroundings. Usually, we control the **Temperature** (using a heat bath) and the **Pressure** (by exposing the system to the atmosphere). 

To simulate these conditions, we change the integrator in LAMMPS.

### The NVT Ensemble (Constant Volume and Temperature)
In the NVT (Canonical) ensemble, the volume of the box is fixed, but the system exchanges heat with an imaginary thermostat to maintain a target temperature. LAMMPS typically uses the **Nosé-Hoover thermostat** for this.



**How to change your script:**
To run an NVT simulation at 300 K, find the `fix` command in Section 4 of your `in.argon` script and replace it with:

```lammps
# Apply the NVT integrator to all atoms
# temp <start_T> <stop_T> <T_damp>
# T_damp determines how rapidly the temperature is relaxed (usually 100x the timestep)
fix             1 all nvt temp 300.0 300.0 100.0

```

*Note:* When using a thermostat, you will notice that your *Total Energy* is no longer conserved—it will fluctuate as the thermostat adds or removes kinetic energy to maintain 300 K.

### The NPT Ensemble (Constant Pressure and Temperature)

In the NPT (Isothermal-Isobaric) ensemble, both the temperature and the pressure are controlled. Because the pressure is held constant, the **volume of the simulation box will expand or contract** dynamically. This is the ensemble you must use if you want to measure the equilibrium density of a liquid or gas at atmospheric conditions.

**How to change your script:**
To run an NPT simulation at 300 K and 1 Atmosphere of pressure, replace the `fix` command with:

```lammps
# Apply the NPT integrator to all atoms
# iso <start_P> <stop_P> <P_damp> controls uniform pressure in all directions
# P_damp is usually 1000x the timestep
fix             1 all npt temp 300.0 300.0 100.0 iso 1.0 1.0 1000.0

```

### Important Rule for NVT and NPT

When running NVT or NPT, the system needs time to adjust to the thermostat and barostat. You must always run an **equilibration phase** before you start collecting your production data. For example:

```lammps
# 1. Equilibrate the system for 5,000 steps
run             5000

# 2. Reset the average thermodynamic counters to zero
reset_timestep  0

# 3. Run the production phase for 10,000 steps to collect data
run             10000

```


```