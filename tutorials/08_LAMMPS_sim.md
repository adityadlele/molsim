# Tutorial 8: LAMMPS Input Files — Building Your First Simulation

**Objective:** Learn to write LAMMPS input scripts by building a simple Argon gas simulation from scratch.

:::{note} Prerequisites
* Completion of Tutorial 6 (MD Engine concepts)
* Basic understanding of Lennard-Jones potential
* Access to Anvil cluster (Tutorial 2)
:::

---

## 1. What is LAMMPS?

**LAMMPS** (Large-scale Atomic/Molecular Massively Parallel Simulator) is one of the most widely-used molecular dynamics codes in the world. Unlike our Python toy models, LAMMPS is production-grade software designed to:

* Simulate millions of atoms efficiently
* Run in parallel across thousands of CPU cores
* Handle complex force fields (proteins, polymers, metals, etc.)
* Integrate with analysis tools and visualization software

LAMMPS is **command-driven**. You write a text file (the "input script") that tells LAMMPS what to do, step by step. This tutorial will teach you the structure of that file.

---

## 2. The Anatomy of a LAMMPS Input File

A LAMMPS input script is read **sequentially** from top to bottom. It is divided into logical sections (though LAMMPS doesn't enforce this — the sections are for human readability).

**The Standard Structure:**

```text
# ========================================
# SECTION 1: INITIALIZATION
# ========================================
# Define units, boundary conditions, atom style

# ========================================
# SECTION 2: SYSTEM DEFINITION
# ========================================
# Create the simulation box and atoms

# ========================================
# SECTION 3: FORCE FIELD
# ========================================
# Define interactions (Lennard-Jones parameters)

# ========================================
# SECTION 4: SIMULATION SETTINGS
# ========================================
# Set timestep, neighbor lists, fixes

# ========================================
# SECTION 5: OUTPUT
# ========================================
# Specify what data to save

# ========================================
# SECTION 6: EXECUTION
# ========================================
# Run the simulation
```

Each section uses specific LAMMPS commands. We'll build them one at a time.

---

## 3. Section 1: Initialization

The first section sets up the "universe" constants that govern the entire simulation.

### 3.1 Units

LAMMPS supports multiple unit systems. For noble gases and simple fluids, we use **`lj`** (Lennard-Jones reduced units), where:
* Energy is in units of $\epsilon$
* Distance is in units of $\sigma$
* Mass is in units of $m$
* Time is in units of $\tau = \sigma\sqrt{m/\epsilon}$

For Argon, the conversion to real units is:
* $\epsilon = 0.238$ kcal/mol (119.8 K in temperature units)
* $\sigma = 3.405$ Å
* $m = 39.948$ amu

**Command:**
```lammps
units lj
```

### 3.2 Atom Style

The `atom_style` command defines what properties each atom has. For a simple monatomic gas, we use `atomic` (just position, velocity, type).

**Command:**
```lammps
atom_style atomic
```

### 3.3 Boundary Conditions

We want **periodic boundaries** in all three directions (x, y, z) so the gas doesn't escape.

**Command:**
```lammps
boundary p p p
```

### 3.4 Complete Section 1

```lammps
# ========================================
# SECTION 1: INITIALIZATION
# ========================================
units           lj
atom_style      atomic
boundary        p p p
```

---

## 4. Section 2: System Definition

Now we create the simulation box and populate it with atoms.

### 4.1 Defining the Box

We need to specify the size of our cubic simulation box. For a dense Argon system, we'll use a box that gives a reduced density $\rho^* = 0.8$ (liquid-like density).

For 2000 atoms at $\rho^* = 0.8$:
$$\rho^* = \frac{N}{\sigma^3 \times V} = \frac{N}{V} \quad \text{(in reduced units)}$$

$$V = \frac{N}{\rho^*} = \frac{2000}{0.8} = 2500 \, \sigma^3$$

$$L = V^{1/3} = 2500^{1/3} \approx 13.57 \, \sigma$$

**Command:**
```lammps
lattice         fcc 0.8
region          simbox block 0 13.57 0 13.57 0 13.57
create_box      1 simbox
```

**Explanation:**
* `lattice fcc 0.8`: Creates a face-centered cubic lattice template with density 0.8
* `region simbox`: Defines a rectangular region named "simbox"
* `create_box 1 simbox`: Creates the box with space for 1 atom type

### 4.2 Creating Atoms

We'll place atoms on the FCC lattice to avoid overlaps.

**Command:**
```lammps
create_atoms    1 box
```

This creates atoms of type 1 throughout the entire box, using the lattice spacing we defined.

### 4.3 Setting Initial Velocities

We need to give atoms random velocities consistent with 500 K (in reduced units, $T^* = k_B T / \epsilon$).

For Argon at 500 K:
$$T^* = \frac{500 \text{ K}}{119.8 \text{ K}} = 4.17$$

**Command:**
```lammps
velocity        all create 4.17 87287 dist gaussian
```

**Explanation:**
* `all`: Apply to all atoms
* `create 4.17`: Generate velocities for T* = 4.17
* `87287`: Random seed (any integer works)
* `dist gaussian`: Use Gaussian (Maxwell-Boltzmann) distribution

### 4.4 Setting Masses

For reduced units, mass is 1.0.

**Command:**
```lammps
mass            1 1.0
```

### 4.5 Complete Section 2

```lammps
# ========================================
# SECTION 2: SYSTEM DEFINITION
# ========================================
lattice         fcc 0.8
region          simbox block 0 13.57 0 13.57 0 13.57
create_box      1 simbox
create_atoms    1 box
mass            1 1.0
velocity        all create 4.17 87287 dist gaussian
```

---

## 5. Section 3: Force Field

This section defines how atoms interact — the heart of the simulation.

### 5.1 Pair Style

We use the Lennard-Jones 12-6 potential with a cutoff.

**Command:**
```lammps
pair_style      lj/cut 3.0
```

This sets a cutoff distance of $3.0\sigma$ (interactions beyond this are ignored).

### 5.2 Pair Coefficients

For a single atom type interacting with itself via LJ, we specify $\epsilon$ and $\sigma$. In reduced units, both are 1.0.

**Command:**
```lammps
pair_coeff      1 1 1.0 1.0 3.0
```

**Syntax:** `pair_coeff type1 type2 epsilon sigma cutoff`

### 5.3 Neighbor List Settings

LAMMPS builds a "neighbor list" of nearby atoms to speed up force calculations. The `neighbor` command controls this.

**Command:**
```lammps
neighbor        0.3 bin
neigh_modify    every 1 delay 0 check yes
```

**Explanation:**
* `0.3 bin`: Add a 0.3σ "skin" around the cutoff
* `every 1`: Rebuild list every timestep
* `delay 0`: Start building immediately
* `check yes`: Check if rebuild is needed

### 5.4 Complete Section 3

```lammps
# ========================================
# SECTION 3: FORCE FIELD
# ========================================
pair_style      lj/cut 3.0
pair_coeff      1 1 1.0 1.0 3.0
neighbor        0.3 bin
neigh_modify    every 1 delay 0 check yes
```

---

## 6. Section 4: Simulation Settings

### 6.1 Timestep

For LJ systems, a typical timestep is $\Delta t^* = 0.005$ (in reduced units).

**Command:**
```lammps
timestep        0.005
```

### 6.2 Thermostat (NVT Ensemble)

We'll use the Nosé-Hoover thermostat to maintain constant temperature.

**Command:**
```lammps
fix             1 all nvt temp 4.17 4.17 0.5
```

**Explanation:**
* `fix 1`: Assign ID "1" to this fix
* `all`: Apply to all atoms
* `nvt`: Nosé-Hoover thermostat
* `temp 4.17 4.17`: Start and end temperature (T* = 4.17)
* `0.5`: Damping parameter (controls coupling strength, in time units)

### 6.3 Complete Section 4

```lammps
# ========================================
# SECTION 4: SIMULATION SETTINGS
# ========================================
timestep        0.005
fix             1 all nvt temp 4.17 4.17 0.5
```

---

## 7. Section 5: Output

We need to save thermodynamic data and optionally trajectory snapshots.

### 7.1 Thermo Output

The `thermo` command controls how often thermodynamic data is printed.

**Command:**
```lammps
thermo_style    custom step temp press pe ke etotal density
thermo          100
```

**Explanation:**
* `thermo_style custom`: Use custom output format
* `step temp press pe ke etotal density`: Columns to print
* `thermo 100`: Print every 100 steps

### 7.2 Trajectory Output (Optional)

To visualize atoms moving, save coordinates periodically.

**Command:**
```lammps
dump            1 all custom 500 argon.lammpstrj id type x y z vx vy vz
dump_modify     1 sort id
```

**Explanation:**
* `dump 1 all custom`: Create dump with ID "1", include all atoms
* `500`: Save every 500 steps
* `argon.lammpstrj`: Output filename
* `id type x y z vx vy vz`: Data to save
* `dump_modify 1 sort id`: Sort atoms by ID for consistent visualization

### 7.3 Redirecting Thermo to File

To save thermo data to a file for analysis:

**Command:**
```lammps
log             thermo.out
```

This writes all thermo output to `thermo.out` instead of just the screen.

### 7.4 Complete Section 5

```lammps
# ========================================
# SECTION 5: OUTPUT
# ========================================
thermo_style    custom step temp press pe ke etotal density
thermo          100
dump            1 all custom 500 argon.lammpstrj id type x y z vx vy vz
dump_modify     1 sort id
log             thermo.out
```

---

## 8. Section 6: Execution

Finally, we tell LAMMPS to actually run the simulation.

### 8.1 Running the Simulation

**Command:**
```lammps
run             50000
```

This runs 50,000 timesteps. With $\Delta t^* = 0.005$, total simulation time is $50000 \times 0.005 = 250$ reduced time units.

### 8.2 Clean Exit

**Command:**
```lammps
write_data      final_state.data
```

This saves the final positions and velocities to a file, which can be used to restart or continue the simulation later.

### 8.3 Complete Section 6

```lammps
# ========================================
# SECTION 6: EXECUTION
# ========================================
run             50000
write_data      final_state.data
```

---

## 9. Complete Input File

Here is the full, working LAMMPS input script:

```lammps
# ========================================
# LAMMPS Input Script: Argon Gas
# Tutorial 7 - Practical Molecular Simulations
# 2000 atoms, T* = 4.17 (500 K), rho* = 0.8
# ========================================

# ========================================
# SECTION 1: INITIALIZATION
# ========================================
units           lj
atom_style      atomic
boundary        p p p

# ========================================
# SECTION 2: SYSTEM DEFINITION
# ========================================
lattice         fcc 0.8
region          simbox block 0 13.57 0 13.57 0 13.57
create_box      1 simbox
create_atoms    1 box
mass            1 1.0
velocity        all create 4.17 87287 dist gaussian

# ========================================
# SECTION 3: FORCE FIELD
# ========================================
pair_style      lj/cut 3.0
pair_coeff      1 1 1.0 1.0 3.0
neighbor        0.3 bin
neigh_modify    every 1 delay 0 check yes

# ========================================
# SECTION 4: SIMULATION SETTINGS
# ========================================
timestep        0.005
fix             1 all nvt temp 4.17 4.17 0.5

# ========================================
# SECTION 5: OUTPUT
# ========================================
thermo_style    custom step temp press pe ke etotal density
thermo          100
dump            1 all custom 500 argon.lammpstrj id type x y z vx vy vz
dump_modify     1 sort id
log             thermo.out

# ========================================
# SECTION 6: EXECUTION
# ========================================
run             50000
write_data      final_state.data
```

---

## 10. Testing Your Input File on Anvil

Before running a full production job, you can test your input file interactively.

### 10.1 Loading LAMMPS Module

LAMMPS on Anvil requires specific modules to be loaded:

```bash
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310
```

### 10.2 Interactive Testing

Log in to Anvil and request an interactive session:

```bash
srun -p shared -A chm250117 --nodes=1 --ntasks=4 --time=00:30:00 --pty /bin/bash
```

Once your interactive job starts:

```bash
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310
lmp -in argon.in
```

This will run LAMMPS and create the output files in your current directory.

### 10.3 Class Examples Directory

Sample LAMMPS input files are available in the shared class directory:

```bash
ls /anvil/projects/x-chm250117/class_examples/
```

You can copy these examples to your workspace:

```bash
cp /anvil/projects/x-chm250117/class_examples/argon_example.in $SCRATCH/
```

---

## 11. Understanding the Output

After running, you'll have several files:

| File | Contents |
|------|----------|
| `thermo.out` | Thermodynamic data (T, P, E) vs. time |
| `argon.lammpstrj` | Trajectory snapshots for visualization |
| `final_state.data` | Final configuration (for restart) |
| `log.lammps` | Full LAMMPS log (includes warnings, timings) |

### 11.1 Quick Check

Open `thermo.out` and look at the columns:

```text
Step Temp Press PotEng KinEng TotEng Density 
0 4.1700000 16.503422 -2.6604462 6.2487500 3.5883038 0.8000000 
100 4.0523344 15.341847 -2.5327691 6.0722148 3.5394457 0.8000000 
200 4.1134589 15.982564 -2.5932117 6.1637866 3.5705749 0.8000000 
...
```

**What to look for:**
* **Temperature:** Should fluctuate around 4.17 (our target)
* **Total Energy:** Should be roughly constant (energy conservation)
* **Density:** Should remain 0.8 (constant volume)

---

## 12. Common LAMMPS Commands Reference

Here's a quick reference for commands you'll use frequently:

| Command | Purpose | Example |
|---------|---------|---------|
| `units` | Set unit system | `units lj` |
| `atom_style` | Define atom properties | `atom_style atomic` |
| `boundary` | Set boundary conditions | `boundary p p p` |
| `lattice` | Define crystal structure | `lattice fcc 0.8` |
| `region` | Define geometric region | `region box block 0 10 0 10 0 10` |
| `create_box` | Create simulation box | `create_box 1 box` |
| `create_atoms` | Add atoms | `create_atoms 1 box` |
| `mass` | Set atom mass | `mass 1 1.0` |
| `velocity` | Set initial velocities | `velocity all create 1.0 12345` |
| `pair_style` | Choose interaction model | `pair_style lj/cut 2.5` |
| `pair_coeff` | Set interaction parameters | `pair_coeff 1 1 1.0 1.0` |
| `timestep` | Set integration timestep | `timestep 0.005` |
| `fix` | Apply constraints/thermostats | `fix 1 all nvt temp 1.0 1.0 0.1` |
| `thermo` | Output frequency | `thermo 100` |
| `dump` | Save trajectory | `dump 1 all custom 1000 traj.lammpstrj` |
| `run` | Execute simulation | `run 10000` |

---

## 13. Exercise: Modify the Temperature

**Task:** Create a second input file that simulates Argon at **300 K** instead of 500 K.

**Steps:**
1. Copy `argon.in` to `argon_300K.in`
2. Calculate the reduced temperature: $T^* = 300 / 119.8 = 2.50$
3. Change the `velocity` and `fix nvt` commands to use 2.50
4. Change the output filenames to avoid overwriting
5. Run both simulations and compare average temperatures

<details>
<summary>Click to see the solution</summary>

**Modified lines:**
```lammps
velocity        all create 2.50 87287 dist gaussian
fix             1 all nvt temp 2.50 2.50 0.5
log             thermo_300K.out
dump            1 all custom 500 argon_300K.lammpstrj id type x y z vx vy vz
```
</details>

---

## 14. Next Steps

Now that you understand the structure of a LAMMPS input file:
* **Tutorial 8:** Running LAMMPS on the Anvil cluster
* **Tutorial 9:** Analyzing LAMMPS output with Python

:::{tip} Best Practices
* **Always comment your input files** — explain what each section does
* **Use descriptive variable names** for regions and fixes
* **Start with short test runs** (1000 steps) before committing to long production runs
* **Check energy conservation** before trusting results
* **Keep a library of working input files** as templates for future projects
:::

---

## 15. Further Reading

* **[LAMMPS Manual](https://docs.lammps.org/)** — The official documentation (searchable!)
* **[Anvil LAMMPS Documentation](https://www.rcac.purdue.edu/knowledge/anvil/software/installing_applications/lammps/provided_module)** — Anvil-specific LAMMPS setup and modules
* **[LAMMPS Tutorial: Argon](https://icme.hpc.msstate.edu/mediawiki/index.php/LAMMPS_Argon)** — Mississippi State tutorial
* **[Molecular Dynamics with LAMMPS](https://lammpstutorials.github.io/)** — Community tutorials
* **Frenkel & Smit, Chapter 4** — Reduced units and their conversions

:::{warning} Units Matter!
The most common mistake in LAMMPS is mixing unit systems. If you use `units real` (real-world units like Angstroms and kcal/mol), you **cannot** use reduced LJ parameters (epsilon=1, sigma=1). Always double-check your units are consistent throughout the input file.
:::