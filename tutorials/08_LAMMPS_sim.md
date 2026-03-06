# Tutorial 8: LAMMPS Input Files — Building Your First Simulation

**Objective:** Learn to write LAMMPS input scripts by building a simple Argon gas simulation from scratch.

:::{note} Prerequisites
* Completion of Tutorial 6 (MD Engine concepts)
* Basic understanding of Lennard-Jones potential
* Access to Anvil cluster (Tutorial 2)
:::

:::{important} In-Class Lab Activity
This tutorial is designed as a hands-on lab. Follow along and execute each step on Anvil. Raise your hand if you get stuck!
:::

---

## 1. Getting Started: Creating Your Workspace

Before we write any LAMMPS code, let's set up a proper workspace on Anvil.

### Step 1: Log in to Anvil

Open a terminal on Anvil:
* **Via SSH:** Use your terminal application
* **Via OnDemand:** Go to [ondemand.anvil.rcac.purdue.edu](https://ondemand.anvil.rcac.purdue.edu) → **Clusters** → **Anvil Shell Access**

### Step 2: Navigate to Your Scratch Space

```bash
cd $SCRATCH
```

:::{tip} Why Scratch?
`$SCRATCH` is a high-speed storage area designed for running simulations. It's much faster than your home directory.
:::

### Step 3: Create a Tutorial Directory

```bash
mkdir lammps_tutorial
cd lammps_tutorial
pwd
```

The `pwd` command shows your current location. You should see:
```text
/anvil/scratch/x-username/lammps_tutorial
```

### Step 4: Create Your First LAMMPS Input File

We'll use the `nano` text editor (simple and beginner-friendly):

```bash
nano argon.in
```

This opens a blank file. You'll see a text editor with a menu at the bottom.

:::{note} Nano Quick Reference
* **Type:** Just start typing
* **Save:** Press `Ctrl+O`, then `Enter`
* **Exit:** Press `Ctrl+X`
* **Cut line:** `Ctrl+K`
* **Paste:** `Ctrl+U`
:::

**Keep nano open** — we'll fill it with LAMMPS commands in the next section!

---

## 2. What is LAMMPS?

**LAMMPS** (Large-scale Atomic/Molecular Massively Parallel Simulator) is one of the most widely-used molecular dynamics codes in the world. Unlike our Python toy models, LAMMPS is production-grade software designed to:

* Simulate millions of atoms efficiently
* Run in parallel across thousands of CPU cores
* Handle complex force fields (proteins, polymers, metals, etc.)
* Integrate with analysis tools and visualization software

LAMMPS is **command-driven**. You write a text file (the "input script") that tells LAMMPS what to do, step by step. This tutorial will teach you the structure of that file.

---

## 2. What is LAMMPS?

**LAMMPS** (Large-scale Atomic/Molecular Massively Parallel Simulator) is production-grade MD software used worldwide. Unlike our Python toy model, LAMMPS can:
* Simulate millions of atoms efficiently
* Run in parallel across thousands of cores
* Handle complex materials (proteins, polymers, metals)

**How it works:** You write a text file (input script) telling LAMMPS what to do, step by step.

---

## 3. Building Your Input File: Section by Section

Now let's fill in the `argon.in` file you have open in nano. Copy each section below.

### The Six-Section Structure

Every LAMMPS input follows this template:

```text
# SECTION 1: INITIALIZATION - Define the "universe"
# SECTION 2: SYSTEM DEFINITION - Create atoms
# SECTION 3: FORCE FIELD - Define interactions
# SECTION 4: SIMULATION SETTINGS - Timestep, thermostat
# SECTION 5: OUTPUT - What to save
# SECTION 6: EXECUTION - Run!
```

---

## 4. Section 1: Initialization

**Copy this into nano:**

```lammps
# ========================================
# SECTION 1: INITIALIZATION
# ========================================
units           real
atom_style      atomic
boundary        p p p
```

**What it does:**
* `units real` — Use real-world units: Angstroms (Å), femtoseconds (fs), kcal/mol, Kelvin (K)
* `atom_style atomic` — Simple atoms (no bonds or molecules)
* `boundary p p p` — Periodic boundaries (atoms wrap around edges)

:::{note} Why Real Units?
Real units let us use actual temperatures (300 K instead of T*=2.50) and directly compare with experimental data. All parameters (distances, energies) are in familiar scientific units.
:::

---

## 5. Section 2: System Definition

**Add this to your file:**

```lammps
# ========================================
# SECTION 2: SYSTEM DEFINITION
# ========================================
# FCC lattice with 3.8 Angstrom conventional cell
lattice         fcc 3.8
# Simulation box: 8 lattice units per edge
region          simbox block 0 8 0 8 0 8
create_box      1 simbox
create_atoms    1 region simbox
mass            1 39.948
velocity        all create 500.0 87287 dist gaussian
```

**What it does:**
* Creates a simulation box with ~2048 Argon atoms
* `lattice fcc 3.8` — FCC lattice with conventional cell size 3.8 Å
* `region simbox block 0 8 0 8 0 8` — Define a region extending 8 lattice units  
* `create_box 1 simbox` — Create simulation box with 1 atom type
* `create_atoms 1 region simbox` — Fill region with atoms on FCC lattice
* `mass 1 39.948` — Argon atomic mass (amu)  
* Gives atoms random velocities for T = 500 K

:::{important} CRITICAL: Lattice Units vs Angstroms
**This is the most confusing part of LAMMPS!**

After using `lattice`, the `region` command uses **lattice units**, NOT Angstroms!

- `lattice fcc 3.8` sets the FCC cell size to 3.8 Å
- `region block 0 8 0 8 0 8` means 8 **lattice units** (= 8 FCC cells)
- **Actual box size in Angstroms:** 8 × 3.8 = **30.4 Å**
- **Total atoms:** 8³ unit cells × 4 atoms/cell = **2048 atoms**

To use Angstroms directly in region, add `units box` before region:
```lammps
lattice         fcc 3.8
region          simbox block 0 30.4 0 30.4 0 30.4 units box
create_box      1 simbox
create_atoms    1 region simbox
```
Both give the same result!
:::

:::{tip} Quick Reference
For ~2000 atoms with FCC:
```lammps
lattice         fcc 3.8
region          simbox block 0 8 0 8 0 8
create_box      1 simbox
create_atoms    1 region simbox    # 8 cells = 2048 atoms
```

For fewer atoms:
```lammps
lattice         fcc 3.8  
region          simbox block 0 5 0 5 0 5
create_box      1 simbox
create_atoms    1 region simbox    # 5 cells = 500 atoms
```

For more atoms:
```lammps
lattice         fcc 3.8
region          simbox block 0 10 0 10 0 10
create_box      1 simbox
create_atoms    1 region simbox  # 10 cells = 4000 atoms
```
:::

---

## 6. Section 3: Force Field

**Add this:**

```lammps
# ========================================
# SECTION 3: FORCE FIELD
# ========================================
pair_style      lj/cut 12.0
pair_coeff      1 1 0.238 3.405
neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes
```

**What it does:**
* `lj/cut 12.0` — Lennard-Jones potential with cutoff at 12.0 Å
* `pair_coeff 1 1 0.238 3.405` — For Argon:
  * ε = 0.238 kcal/mol (well depth)
  * σ = 3.405 Å (atom diameter)
* `neighbor 2.0` — Build neighbor list with 2.0 Å buffer

:::{note} Argon LJ Parameters
These are the standard Lennard-Jones parameters for Argon from experimental fitting:
* ε = 0.238 kcal/mol = 119.8 K × k_B
* σ = 3.405 Å
These reproduce experimental properties of Argon accurately.
:::

---

## 7. Section 4: Simulation Settings

**Add this:**

```lammps
# ========================================
# SECTION 4: SIMULATION SETTINGS
# ========================================
timestep        2.0
fix             1 all nvt temp 500.0 500.0 100.0
```

**What it does:**
* `timestep 2.0` — Time advances in steps of 2.0 fs (femtoseconds)
* `fix 1 all nvt temp 500.0 500.0 100.0` — Nosé-Hoover thermostat at 500 K
  * Start temp: 500 K
  * End temp: 500 K
  * Damping: 100 fs

:::{note} Timestep Selection
2.0 fs is standard for simple Lennard-Jones systems. For molecules with bonds, you'd use 1.0 fs.
:::

---

## 8. Section 5: Output

**Add this:**

```lammps
# ========================================
# SECTION 5: OUTPUT
# ========================================
thermo_style    custom step temp press pe ke etotal density
thermo          500
dump            1 all custom 1000 argon.lammpstrj id type x y z vx vy vz
dump_modify     1 sort id
log             thermo.out
```

**What it does:**
* `thermo 500` — Print thermodynamic data every 500 steps (= 1 ps)
* `dump` — Save atom positions/velocities every 1000 steps (= 2 ps)
* `log thermo.out` — Write output to file (for later analysis)

---

## 9. Section 6: Execution

**Final section:**

```lammps
# ========================================
# SECTION 6: EXECUTION
# ========================================
run             50000
write_data      final_state.data
```

**What it does:**
* `run 50000` — Simulate for 50,000 timesteps = 100 ps (picoseconds)
* `write_data` — Save final configuration (for restart/continuation)

:::{tip} Simulation Time
50,000 steps × 2 fs/step = 100,000 fs = 100 ps = 0.1 ns

This is enough to equilibrate Argon and get basic statistics.
:::

---

## 10. Save and Exit

Now save your file:
1. Press `Ctrl+O` (that's the letter O, not zero)
2. Press `Enter` to confirm filename
3. Press `Ctrl+X` to exit nano

**Verify your file was created:**
```bash
ls -lh argon.in
```

You should see:
```text
-rw-r--r-- 1 username group 1.2K date argon.in
```

---

## 11. Creating an Organized Folder Structure

Professional MD simulations require good organization. Let's create a proper folder structure for running multiple simulations.

### Why Folder Structure Matters

When you run multiple simulations (different temperatures, ensembles, parameters), you need:
* ✅ Easy to find specific results
* ✅ No accidentally overwriting files
* ✅ Clear documentation of what each simulation does

### The Standard Structure

We'll create this hierarchy:
```text
lammps_tutorial/
├── NVT_simulations/
│   ├── T_300K/
│   ├── T_500K/
│   └── T_700K/
└── NPT_simulations/
    ├── P_1atm/
    └── P_10atm/
```

### Step 1: Create the Main Directories

```bash
cd $SCRATCH/lammps_tutorial

# Create top-level folders
mkdir -p NVT_simulations
mkdir -p NPT_simulations

# Verify
ls -lh
```

You should see:
```text
drwxr-xr-x  NVT_simulations/
drwxr-xr-x  NPT_simulations/
-rw-r--r--  argon.in
```

### Step 2: Create Temperature Subfolders (NVT)

```bash
cd NVT_simulations

# Create three temperature folders
mkdir T_300K T_500K T_700K

# Check
ls
```

Output:
```text
T_300K  T_500K  T_700K
```

### Step 3: Copy Input Files for Each Temperature

We'll create three versions of the input file, one per temperature.

**For T_300K:**
```bash
cd T_300K
nano argon_300K.in
```

Paste this **complete** input file:

```lammps
# ========================================
# Argon NVT Simulation at 300 K
# ========================================

# ========================================
# SECTION 1: INITIALIZATION
# ========================================
units           real
atom_style      atomic
boundary        p p p

# ========================================
# SECTION 2: SYSTEM DEFINITION
# ========================================
lattice         fcc 3.8
region          simbox block 0 8 0 8 0 8
create_box      1 simbox
create_atoms    1 region simbox
mass            1 39.948
velocity        all create 300.0 87287 dist gaussian

# ========================================
# SECTION 3: FORCE FIELD
# ========================================
pair_style      lj/cut 12.0
pair_coeff      1 1 0.238 3.405
neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes

# ========================================
# SECTION 4: SIMULATION SETTINGS
# ========================================
timestep        2.0
fix             1 all nvt temp 300.0 300.0 100.0

# ========================================
# SECTION 5: OUTPUT
# ========================================
thermo_style    custom step temp press pe ke etotal density
thermo          500
dump            1 all custom 1000 argon_300K.lammpstrj id type x y z vx vy vz
dump_modify     1 sort id
log             thermo_300K.out

# ========================================
# SECTION 6: EXECUTION
# ========================================
run             50000
write_data      final_300K.data
```

Save: `Ctrl+O`, `Enter`, `Ctrl+X`

**For T_500K:**
```bash
cd ../T_500K
nano argon_500K.in
```

Paste the same file BUT change these lines:
```lammps
# Header comment: "Argon NVT Simulation at 500 K"
velocity        all create 500.0 87287 dist gaussian
fix             1 all nvt temp 500.0 500.0 100.0
dump            1 all custom 1000 argon_500K.lammpstrj id type x y z vx vy vz
log             thermo_500K.out
write_data      final_500K.data
```

**For T_700K:**
```bash
cd ../T_700K
nano argon_700K.in
```

Change to 700 K:
```lammps
# Header comment: "Argon NVT Simulation at 700 K"
velocity        all create 700.0 87287 dist gaussian
fix             1 all nvt temp 700.0 700.0 100.0
dump            1 all custom 1000 argon_700K.lammpstrj id type x y z vx vy vz
log             thermo_700K.out
write_data      final_700K.data
```

### Step 4: Verify NVT Folder Structure

```bash
cd $SCRATCH/lammps_tutorial/NVT_simulations
tree
```

Or without tree:
```bash
ls -R
```

You should see:
```text
T_300K/:
argon_300K.in

T_500K/:
argon_500K.in

T_700K/:
argon_700K.in
```

---

## 12. Creating NPT Simulation Inputs

NPT simulations control both temperature AND pressure. Let's create two pressure conditions.

### Step 1: Create Pressure Subfolders

```bash
cd $SCRATCH/lammps_tutorial/NPT_simulations

mkdir P_1atm P_10atm

ls
```

Output:
```text
P_1atm  P_10atm
```

### Step 2: NPT Input for 1 atm

```bash
cd P_1atm
nano argon_npt_1atm.in
```

Paste this:

```lammps
# ========================================
# Argon NPT Simulation at 1 atm
# T = 500 K, P = 1 atm
# ========================================

# ========================================
# INITIALIZATION
# ========================================
units           real
atom_style      atomic
boundary        p p p

# ========================================
# SYSTEM DEFINITION
# ========================================
lattice         fcc 3.8
region          simbox block 0 8 0 8 0 8
create_box      1 simbox
create_atoms    1 region simbox
mass            1 39.948
velocity        all create 500.0 87287 dist gaussian

# ========================================
# FORCE FIELD
# ========================================
pair_style      lj/cut 12.0
pair_coeff      1 1 0.238 3.405
neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes

# ========================================
# SIMULATION SETTINGS (NPT!)
# ========================================
timestep        2.0
fix             1 all npt temp 500.0 500.0 100.0 iso 1.0 1.0 1000.0

# ========================================
# OUTPUT
# ========================================
thermo_style    custom step temp press pe ke etotal vol density
thermo          500
dump            1 all custom 2000 argon_npt_1atm.lammpstrj id type x y z
dump_modify     1 sort id
log             thermo_npt_1atm.out

# ========================================
# EXECUTION
# ========================================
run             100000
write_data      final_npt_1atm.data
```

:::{note} NPT Key Differences
* `fix npt` instead of `fix nvt`
* `iso 1.0 1.0 1000.0` — isotropic pressure control at 1 atm
  * Start pressure: 1.0 atm
  * End pressure: 1.0 atm
  * Damping: 1000.0 fs (pressure control is slower than temperature!)
* `vol` added to thermo output (box volume changes!)
* Longer run (100,000 steps = 200 ps) — pressure equilibration is slower
:::

:::{tip} Pressure Units
In `units real`, pressure is in atmospheres (atm). 1 atm = 101.325 kPa = standard atmospheric pressure at sea level.
:::

:::{important} Initial Box Size for NPT
NPT simulations start with the same box size as NVT (21.04 Å), but the box will shrink or expand during the simulation to maintain the target pressure. At 500 K and 1 atm, Argon is a gas, so expect the box to expand significantly!
:::

### Step 3: NPT Input for 10 atm

```bash
cd ../P_10atm
nano argon_npt_10atm.in
```

Copy the 1 atm file but change:
```lammps
# Header: "Argon NPT Simulation at 10 atm"
# T = 500 K, P = 10 atm
fix             1 all npt temp 500.0 500.0 100.0 iso 10.0 10.0 1000.0
dump            1 all custom 2000 argon_npt_10atm.lammpstrj id type x y z
log             thermo_npt_10atm.out
write_data      final_npt_10atm.data
```

### Step 4: Verify Complete Structure

```bash
cd $SCRATCH/lammps_tutorial
tree
```

Or:
```bash
find . -name "*.in" -type f
```

You should see all 5 input files:
```text
./NVT_simulations/T_300K/argon_300K.in
./NVT_simulations/T_500K/argon_500K.in
./NVT_simulations/T_700K/argon_700K.in
./NPT_simulations/P_1atm/argon_npt_1atm.in
./NPT_simulations/P_10atm/argon_npt_10atm.in
```

---

## 13. Understanding Real Units vs Reduced Units

### Why We Use Real Units in This Class

**Real units** (`units real` in LAMMPS):
- ✅ Direct: Temperature = 300 K (not T* = 2.50)
- ✅ Intuitive: Pressure = 1 atm (not P* = 0.002)
- ✅ Comparable: Can directly compare with experimental data
- ✅ Practical: Easy to understand physical meaning

**Reduced units** (`units lj`):
- Good for theory and dimensionless comparisons
- Common in older MD papers
- Less intuitive for beginners

### Units in LAMMPS Real Mode

| Quantity | Units | Example |
|----------|-------|---------|
| Distance | Angstroms (Å) | σ = 3.405 Å |
| Time | femtoseconds (fs) | timestep = 2.0 fs |
| Energy | kcal/mol | ε = 0.238 kcal/mol |
| Temperature | Kelvin (K) | T = 500 K |
| Pressure | atmospheres (atm) | P = 1 atm |
| Mass | atomic mass units (amu) | m_Ar = 39.948 amu |

### Argon Parameters (Memorize These!)

| Parameter | Value | Physical Meaning |
|-----------|-------|------------------|
| ε | 0.238 kcal/mol | LJ well depth |
| σ | 3.405 Å | LJ atom diameter |
| mass | 39.948 amu | Atomic mass |
| Cutoff | 12.0 Å | ~3.5σ (interaction range) |

:::{tip} Converting to Reduced Units (Optional)
If you see papers using reduced units:
* T* = T / T_ε, where T_ε = ε/k_B = 119.8 K
* So 300 K = 300/119.8 = 2.50 (reduced)
* So 500 K = 500/119.8 = 4.17 (reduced)

But you don't need to do this in our class!
:::

---

## 14. Summary Checklist

After completing this tutorial, you should have:

- [ ] Created organized folder structure:
  - `NVT_simulations/` with 3 temperature folders
  - `NPT_simulations/` with 2 pressure folders
- [ ] 5 complete LAMMPS input files (.in files)
- [ ] Understand the difference between NVT and NPT fixes
- [ ] Know how to convert real units (K, atm) to reduced units

**File count check:**
```bash
cd $SCRATCH/lammps_tutorial
find . -name "*.in" | wc -l
```

Should output: `5`

**Next Tutorial:** Running these simulations on Anvil (Tutorial 8)

---

## 15. Troubleshooting

**Problem: "mkdir: cannot create directory: File exists"**
* That folder already exists. Use `ls` to check what's there, or use `-p` flag: `mkdir -p foldername`

**Problem: "No such file or directory"**
* You're in the wrong directory. Use `pwd` to check location, `cd` to navigate

**Problem: Lost track of where you are**
```bash
cd $SCRATCH/lammps_tutorial
pwd
ls -R
```

**Problem: Want to start over**
```bash
cd $SCRATCH
rm -rf lammps_tutorial
# Then start from Section 11, Step 1
```

:::{warning} Be Careful with `rm -rf`
This permanently deletes everything! Double-check you're in the right directory first with `pwd`.
:::

---

## 16. Class Examples

Pre-made complete folder structures are available:
```bash
ls /anvil/projects/x-chm250117/class_examples/
```

If you get stuck, you can copy the complete structure:
```bash
cd $SCRATCH
cp -r /anvil/projects/x-chm250117/class_examples/lammps_tutorial ./
```

---

## 17. Further Reading

* **[LAMMPS Manual](https://docs.lammps.org/)** — Command reference
* **[Anvil LAMMPS Documentation](https://www.rcac.purdue.edu/knowledge/anvil/software/installing_applications/lammps/provided_module)** — Anvil-specific guide
* **[LAMMPS fix npt](https://docs.lammps.org/fix_nh.html)** — NPT thermostat/barostat documentation

Let's make sure your input file works before submitting a job.

### Step 1: Request Interactive Session

```bash
srun -p shared -A chm250117 --nodes=1 --ntasks=4 --time=00:10:00 --pty /bin/bash
```

Wait for the prompt to change (30 seconds to 2 minutes). When you see a new prompt, you're on a compute node.

### Step 2: Load LAMMPS

```bash
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310
```

### Step 3: Run LAMMPS

```bash
lmp -in argon.in
```

You should see output scrolling by:
```text
LAMMPS (10 Mar 2021)
...
Step Temp Press PotEng KinEng TotEng Density
0 4.1700000 16.503422 -2.6604462 6.2487500 3.5883038 0.8000000
100 4.0523344 15.341847 -2.5327691 6.0722148 3.5394457 0.8000000
...
```

**If you see this, success!** Your input file works.

### Step 4: Exit Interactive Session

```bash
exit
```

---

## 12. Understanding the Output Files

## 12. Understanding the Output Files

After running, check what files were created:

```bash
ls -lh
```

| File | Contents |
|------|----------|
| `thermo.out` | Temperature, pressure, energy vs. time |
| `argon.lammpstrj` | Trajectory snapshots (for visualization) |
| `final_state.data` | Final positions (for restart) |
| `log.lammps` | Full LAMMPS log (warnings, timings) |

**Quick check:**
```bash
head -20 thermo.out
```

Look for the temperature column — it should fluctuate around 4.17 (your target).

---

## 13. Common LAMMPS Commands Reference

| Command | Purpose | Example |
|---------|---------|---------|
| `units` | Set unit system | `units lj` |
| `boundary` | Set boundary conditions | `boundary p p p` |
| `create_box` | Create simulation box | `create_box 1 simbox` |
| `create_atoms` | Add atoms | `create_atoms 1 box` |
| `mass` | Set atom mass | `mass 1 1.0` |
| `velocity` | Set initial velocities | `velocity all create 1.0 12345` |
| `pair_style` | Choose interaction model | `pair_style lj/cut 2.5` |
| `pair_coeff` | Set interaction parameters | `pair_coeff 1 1 1.0 1.0` |
| `timestep` | Set integration timestep | `timestep 0.005` |
| `fix` | Apply thermostat/barostat | `fix 1 all nvt temp 1.0 1.0 0.1` |
| `thermo` | Output frequency | `thermo 100` |
| `dump` | Save trajectory | `dump 1 all custom 1000 traj.lammpstrj` |
| `run` | Execute simulation | `run 10000` |

---

## 14. Exercise: Modify the Temperature

**Challenge:** Create a second input file for Argon at **300 K** instead of 500 K.

<details>
<summary>Hint 1: What needs to change?</summary>

Calculate the reduced temperature:
$$T^* = \frac{T}{T_\epsilon} = \frac{300 \text{ K}}{119.8 \text{ K}} = 2.50$$

You need to change TWO lines in the input file.
</details>

<details>
<summary>Hint 2: Which lines?</summary>

1. The `velocity` line (sets initial temperature)
2. The `fix nvt` line (thermostat target)
</details>

<details>
<summary>Solution</summary>

```bash
# Copy the file
cp argon.in argon_300K.in

# Edit it
nano argon_300K.in
```

Change these two lines:
```lammps
velocity        all create 2.50 87287 dist gaussian
fix             1 all nvt temp 2.50 2.50 0.5
```

Also update the output filenames to avoid overwriting:
```lammps
dump            1 all custom 500 argon_300K.lammpstrj id type x y z vx vy vz
log             thermo_300K.out
```

Test it in an interactive session!
</details>

---

## 15. Class Examples Directory

Pre-made example files are available:

```bash
ls /anvil/projects/x-chm250117/class_examples/
```

Copy an example to your directory:
```bash
cp /anvil/projects/x-chm250117/class_examples/argon_example.in ./
```

---

## 16. Troubleshooting

**Problem: "nano: command not found"**
* You're not on Anvil or in a proper terminal

**Problem: "Permission denied"**
* Make sure you're in `$SCRATCH`, not `$HOME`

**Problem: LAMMPS crashes with "Unknown command"**
* Check for typos in command names
* Make sure there are no extra spaces at line beginnings

**Problem: "ERROR: Cannot open input script"**
* Check filename: `ls -lh argon.in`
* Make sure you're in the right directory: `pwd`

---

## 17. Summary Checklist

After completing this tutorial, you should have:

- [ ] Created `argon.in` in `$SCRATCH/lammps_tutorial/`
- [ ] Successfully run LAMMPS in an interactive session
- [ ] Seen `thermo.out` and `argon.lammpstrj` files created
- [ ] Understand the 6-section structure of LAMMPS input files

**Next Tutorial:** Running LAMMPS jobs on the cluster (Tutorial 8)

---

## 18. Further Reading

* **[LAMMPS Manual](https://docs.lammps.org/)** — Searchable command reference
* **[Anvil LAMMPS Documentation](https://www.rcac.purdue.edu/knowledge/anvil/software/installing_applications/lammps/provided_module)** — Anvil-specific setup
* **[LAMMPS Tutorials](https://lammpstutorials.github.io/)** — Community tutorials

:::{warning} Units Matter!
The most common mistake is mixing unit systems. If you use `units real`, you **cannot** use ε=1, σ=1. Always check that all parameters match your chosen unit system.
:::. It is divided into logical sections (though LAMMPS doesn't enforce this — the sections are for human readability).

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
create_atoms    1 region simbox
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
create_atoms    1 region simbox
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
create_atoms    1 region simbox
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