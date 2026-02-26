# Practical Molecular Dynamics with LAMMPS

**Objective:** Understand the practical implementation of molecular dynamics simulations, from theory to executable code, using LAMMPS as the implementation framework.

---

## 1. Introduction: From Theory to Practice

In previous lectures, we've covered the theoretical foundations of molecular dynamics:
* Hamilton's equations of motion
* Force fields and interatomic potentials
* Integration algorithms (Velocity Verlet)
* Thermostats and barostats
* Ensemble theory

Now we must translate this theory into actual simulations. This requires understanding:
1. How MD codes are structured
2. What decisions you must make as the researcher
3. How to interpret and validate results

:::{note} The Philosophy
Writing your own MD code (as in Tutorial 6) teaches you *what* is happening inside the black box. Using production MD software (LAMMPS, GROMACS, NAMD) teaches you *how* to apply MD to real research problems. Both skills are essential.
:::

---

## 2. The MD Software Landscape

### 2.1 Major MD Packages

| Software | Strengths | Typical Applications |
|----------|-----------|---------------------|
| **LAMMPS** | Materials, generality, parallelization | Metals, polymers, granular materials, simple fluids |
| **GROMACS** | Biomolecules, speed, analysis tools | Proteins, lipids, drug discovery |
| **NAMD** | Biomolecules, GPU acceleration | Large biomolecular systems, membrane proteins |
| **AMBER** | Force fields for biomolecules | DNA, RNA, proteins |
| **VASP/CP2K** | Ab initio MD (quantum mechanics) | Reactive chemistry, electronic structure |

**Why LAMMPS for this course?**
* **Open source and free** — No licensing barriers
* **Well-documented** — Extensive manual and active user community
* **General-purpose** — Can simulate almost any atomic system
* **HPC-ready** — Designed for parallel supercomputers
* **Extensible** — Easy to add custom modifications

### 2.2 Common Architecture

Despite differences, most MD codes share a common structure:

```text
┌─────────────────────────────────────┐
│   1. INITIALIZATION                 │
│   - Set units, boundary conditions  │
│   - Define atom types               │
└──────────────┬──────────────────────┘
               ↓
┌─────────────────────────────────────┐
│   2. SYSTEM BUILDING                │
│   - Create simulation box           │
│   - Add atoms (from file or lattice)│
│   - Set initial velocities          │
└──────────────┬──────────────────────┘
               ↓
┌─────────────────────────────────────┐
│   3. FORCE FIELD                    │
│   - Define interactions (LJ, bonds) │
│   - Set cutoffs and parameters      │
└──────────────┬──────────────────────┘
               ↓
┌─────────────────────────────────────┐
│   4. SIMULATION PROTOCOL            │
│   - Set timestep                    │
│   - Choose ensemble (NVE/NVT/NPT)   │
│   - Define output frequency         │
└──────────────┬──────────────────────┘
               ↓
┌─────────────────────────────────────┐
│   5. EXECUTION                      │
│   - Energy minimization (optional)  │
│   - Equilibration run               │
│   - Production run                  │
└──────────────┬──────────────────────┘
               ↓
┌─────────────────────────────────────┐
│   6. ANALYSIS                       │
│   - Extract thermodynamic averages  │
│   - Calculate structural properties │
│   - Visualize trajectories          │
└─────────────────────────────────────┘
```

---

## 3. The LAMMPS Input Script Philosophy

### 3.1 Sequential Execution

Unlike compiled programs where you define functions and then call them, LAMMPS executes commands **line-by-line, top to bottom**. Think of it as a recipe:

```text
Step 1: Preheat oven to 350°F    →  units lj
Step 2: Mix flour and sugar      →  create_atoms 1 box
Step 3: Add eggs                 →  velocity all create 1.0 12345
Step 4: Bake for 30 minutes      →  run 50000
```

This means:
* **Order matters** — You can't use a variable before defining it
* **State is persistent** — Settings from earlier commands remain active
* **You can modify on-the-fly** — Can change thermostats mid-simulation

### 3.2 The Command Structure

Every LAMMPS command has a specific syntax:

```lammps
command_name  argument1  argument2  keyword1 value1  keyword2 value2
```

**Example:**
```lammps
fix  1  all  nvt  temp 300.0 300.0  100.0
│    │   │    │    │         │       │
│    │   │    │    └─────────┴───────┴─ Required arguments
│    │   │    └───────────────────────── Fix style (Nosé-Hoover)
│    │   └────────────────────────────── Atom group to apply fix
│    └────────────────────────────────── Fix ID (user-defined name)
└─────────────────────────────────────── Command name
```

### 3.3 Comments and Formatting

```lammps
# This is a comment (ignored by LAMMPS)

units lj                    # Inline comments are allowed

# Blank lines are ignored - use them to organize sections

# Long commands can be split with &
pair_coeff  1 1  1.0 1.0 &
            2.5              # Cutoff on next line
```

---

## 4. Critical Decision Points

### 4.1 Unit Systems

LAMMPS supports multiple unit systems. **This is the #1 source of errors.**

| Unit System | Length | Time | Energy | Mass | Temperature |
|-------------|--------|------|--------|------|-------------|
| **lj** | σ | τ | ε | m | $k_B T/\epsilon$ |
| **real** | Å | fs | kcal/mol | amu | K |
| **metal** | Å | ps | eV | amu | K |
| **si** | m | s | J | kg | K |

**Golden Rule:** All parameters in your input file must be **consistent** with the chosen unit system.

**Example Error:**
```lammps
units real                           # Using real units (kcal/mol)
pair_coeff 1 1  1.0 1.0 2.5          # But these are LJ reduced units!
                                     # ✗ WRONG - Will give nonsense results
```

**Correct:**
```lammps
units real
pair_coeff 1 1  0.238 3.405 10.0     # Argon in real units (kcal/mol, Å)
```

### 4.2 Reduced vs. Real Units

For simple systems (noble gases, generic polymers), **reduced LJ units** are preferred:
* All parameters are dimensionless (ε = 1, σ = 1, m = 1)
* Results are general — apply to any LJ system
* Easy to compare with theory and literature

For specific materials, use **real units**:
* Energy in kcal/mol, distance in Å
* Can directly compare with experimental data
* Required when using biomolecular force fields (AMBER, CHARMM)

**Conversion Example (Argon):**

In reduced units:
```lammps
units lj
pair_coeff 1 1  1.0 1.0 3.0          # ε=1, σ=1, cutoff=3σ
velocity all create 1.0 12345        # T* = 1.0
```

In real units:
```lammps
units real
pair_coeff 1 1  0.238 3.405 10.215   # ε in kcal/mol, σ in Å
velocity all create 119.8 12345      # T = 119.8 K (same as T*=1.0)
```

Where:
* $\epsilon = 0.238$ kcal/mol = 119.8 K × $k_B$
* $\sigma = 3.405$ Å
* $T^* = 1.0$ corresponds to $T = \epsilon/k_B = 119.8$ K

---

## 5. Ensemble Selection

### 5.1 NVE: The Fundamental Ensemble

**Microcanonical Ensemble** — Constant N, V, E

```lammps
timestep 0.005
run 100000                           # No thermostat/barostat = NVE
```

**Uses:**
* Testing energy conservation (numerical stability check)
* Isolated systems (no heat exchange)
* Baseline for comparing thermostats

**What to check:**
* Total energy should fluctuate but not drift
* If energy explodes → timestep too large or bad initial configuration
* If energy drifts slowly → poor integrator or numerical precision issues

### 5.2 NVT: Constant Temperature

**Canonical Ensemble** — Constant N, V, T

```lammps
fix 1 all nvt temp 300.0 300.0 100.0
```

**Parameters:**
* `temp 300.0 300.0` — Start and end temperature (can ramp: `temp 100.0 300.0`)
* `100.0` — Damping parameter $\tau_T$ (in time units)
  * Too small (<10 timesteps): Over-damped, suppresses fluctuations
  * Too large (>1000 timesteps): Under-damped, slow equilibration
  * **Rule of thumb:** $\tau_T \approx 100 \times \Delta t$

**When to use:**
* Most common for production runs
* Studying temperature-dependent properties
* Systems in thermal contact with reservoir (lab experiments)

**Alternatives:**
```lammps
fix 1 all langevin 300.0 300.0 100.0 12345    # Langevin dynamics
fix 1 all temp/berendsen 300.0 300.0 100.0    # Berendsen (equilibration only!)
```

### 5.3 NPT: Constant Pressure

**Isothermal-Isobaric Ensemble** — Constant N, P, T

```lammps
fix 1 all npt temp 300.0 300.0 100.0 iso 1.0 1.0 1000.0
```

**Parameters:**
* `temp 300.0 300.0 100.0` — Temperature control (same as NVT)
* `iso 1.0 1.0 1000.0` — Isotropic pressure control
  * Start pressure: 1.0 atm
  * End pressure: 1.0 atm
  * Damping: 1000.0 (in time units)

**Pressure control modes:**
* `iso` — Uniform scaling (most common for fluids)
* `aniso` — Independent scaling in x, y, z (for crystals under stress)
* `tri` — Full stress tensor (for shear simulations)

**When to use:**
* Determining equilibrium density
* Phase transitions (liquid-vapor, crystal structures)
* Matching experimental conditions (atmospheric pressure)

**Critical:** 
* Pressure control is **much slower** than temperature control
* Requires longer equilibration (10-100× longer than NVT)
* Box size will fluctuate — this is correct!

---

## 6. Timestep Selection

### 6.1 The Stability Criterion

The timestep $\Delta t$ must be **much smaller** than the fastest motion in the system.

**Physical basis:** Atomic vibrations have a period $T_{vib}$. To accurately resolve this motion:
$$\Delta t \ll T_{vib}$$

Typically: $\Delta t \approx T_{vib}/10$

### 6.2 Timestep Guidelines

| System Type | Fastest Motion | Typical $\Delta t$ |
|-------------|----------------|-------------------|
| **Argon (no bonds)** | LJ vibrations | 5 fs (real), 0.005 (LJ) |
| **Water (rigid bonds)** | Constrained H-O | 2 fs |
| **Water (flexible bonds)** | O-H stretch | 0.5 fs |
| **Proteins (SHAKE)** | Constrained bonds | 2 fs |
| **Proteins (flexible)** | C-H stretch | 1 fs |
| **Coarse-grained** | Soft potentials | 10-50 fs |

**In LAMMPS:**
```lammps
# Reduced LJ units
timestep 0.005                       # 0.005 τ

# Real units (femtoseconds)
timestep 1.0                         # 1 fs
```

### 6.3 Testing Your Timestep

**Method:** Run short NVE simulations with different timesteps and plot total energy.

```lammps
# Test run
timestep 0.001
fix 1 all nve
run 10000
write_data test_0.001.data

# Repeat with 0.005, 0.01, 0.02...
```

**Good timestep:** Total energy fluctuates around a constant value (±0.01%)
**Too large:** Total energy drifts upward (exponentially growing errors)

---

## 7. Output and Data Collection

### 7.1 Thermodynamic Output

```lammps
thermo_style custom step temp press pe ke etotal vol density
thermo 100                           # Print every 100 steps
```

**Column selection:**
* **Always include:** `step`, `temp`, `pe`, `ke`, `etotal`
* **For NVT/NPT:** `press`, `vol`, `density`
* **For analysis:** Custom computes (RDF, MSD, etc.)

**Writing to file:**
```lammps
log thermo.out                       # Redirect output to file
```

### 7.2 Trajectory Output (Snapshots)

```lammps
dump 1 all custom 1000 traj.lammpstrj id type x y z vx vy vz
dump_modify 1 sort id                # Sort by atom ID for consistency
```

**Frequency considerations:**
* **Too frequent** (every step): Enormous file sizes, disk I/O bottleneck
* **Too sparse** (every 100,000 steps): Miss important events
* **Rule of thumb:** Every 500-1000 steps for visualization, 5000-10000 for analysis

**File formats:**
* `custom` — LAMMPS native format (text, human-readable)
* `xyz` — Simple format for visualization (VMD, OVITO)
* `dcd` — Binary format (smaller files, faster I/O)

### 7.3 Restart Files

```lammps
restart 50000 restart.*.data         # Save restart every 50k steps
write_restart final.restart          # Save at end
```

**Use cases:**
* Continue a crashed simulation
* Extend a finished simulation
* Branch simulations from the same starting point

**Restarting:**
```lammps
read_restart final.restart
# Continue with new commands
run 100000
```

---

## 8. The Multi-Stage Simulation Protocol

Real research simulations are rarely a single `run` command. The standard workflow has multiple stages:

### 8.1 Stage 1: Energy Minimization

**Purpose:** Remove steric clashes (overlapping atoms)

```lammps
minimize 1.0e-4 1.0e-6 1000 10000
#        etol   ftol   maxiter maxeval
```

**When to use:**
* Always, if building from scratch (lattice or random)
* After modifying a structure (adding molecules, mutations)
* Skip if reading from pre-equilibrated configuration

**Output to check:**
```text
Minimization converged in 547 steps
Final energy: -8342.3 kcal/mol
```

If minimization fails (doesn't converge):
* Initial configuration is very bad (major overlaps)
* Force field parameters are wrong
* Box is too small

### 8.2 Stage 2: NVT Equilibration

**Purpose:** Heat system to target temperature

```lammps
# Start from 0 K (minimized structure has no velocities)
velocity all create 0.0 12345

# Ramp temperature from 0 to 300 K over 50 ps
fix 1 all nvt temp 0.0 300.0 100.0
run 10000                            # 50 ps at dt=0.005 ps
unfix 1                              # Remove this thermostat
```

**Duration:** 10,000-100,000 steps depending on system size

### 8.3 Stage 3: NPT Equilibration (if needed)

**Purpose:** Equilibrate density at target pressure

```lammps
fix 1 all npt temp 300.0 300.0 100.0 iso 1.0 1.0 1000.0
run 100000                           # 500 ps - pressure is slow!
unfix 1
```

**What to monitor:**
* Volume/density should stabilize (no upward/downward trend)
* Pressure oscillates around target (large fluctuations are normal)

### 8.4 Stage 4: Production Run

**Purpose:** Collect data for analysis

```lammps
# Switch to final ensemble (e.g., NVT for canonical averages)
fix 1 all nvt temp 300.0 300.0 100.0

# Reset counters and start fresh output
reset_timestep 0
log production.out

# Long run for statistics
run 1000000                          # 5 ns

write_data final_production.data
```

**Duration:** System-dependent
* Diffusion coefficients: 1-10 ns
* Structural properties: 0.1-1 ns
* Phase transitions: 10-100 ns
* Rare events: microseconds (need enhanced sampling)

---

## 9. Practical Considerations for HPC

### 9.1 Parallel Efficiency

LAMMPS uses **spatial decomposition** — the simulation box is divided among processors.

**Optimal scaling:**
* Each subdomain should have ≥1000 atoms
* For 10,000 atoms: Use 4-16 cores (good efficiency)
* For 100,000 atoms: Use 16-64 cores
* For 1,000,000 atoms: Use 64-256 cores

**Beyond this:** Communication overhead dominates (diminishing returns)

**Test scaling:**
```bash
# Run with different core counts
mpirun -np 4 lmp -in input.in        # Record walltime
mpirun -np 8 lmp -in input.in
mpirun -np 16 lmp -in input.in
# Ideal: 2x cores = 2x speedup (rarely achieved)
```

### 9.2 I/O Optimization

**Disk writes are slow** — they can dominate runtime for large systems.

**Strategies:**
* Write trajectories to fast scratch space (`$SCRATCH`)
* Reduce dump frequency (1000-5000 steps)
* Use binary formats (`dump dcd`) instead of text
* Disable unnecessary thermo output during equilibration

**Example:**
```lammps
# During equilibration - minimal output
thermo 10000
dump 1 all custom 50000 eq.lammpstrj id type x y z

# During production - detailed output
thermo 1000
dump 2 all dcd 1000 prod.dcd
```

### 9.3 Memory Management

**LAMMPS memory scales with:**
* Number of atoms (primary)
* Neighbor list size (cutoff-dependent)
* Number of dumps/computes stored in memory

**Typical usage:** ~100-200 bytes/atom

**For 1 million atoms:** ~200 MB RAM (plus overhead)

**If you run out of memory:**
* Reduce neighbor list skin (`neighbor 0.3 bin` → `neighbor 0.1 bin`)
* Reduce dump frequency
* Don't accumulate time-averaged data over entire run

---

## 10. Common Pitfalls and Debugging

### 10.1 Energy Explosion

**Symptom:** Total energy increases exponentially, simulation crashes

**Causes:**
1. Timestep too large → Reduce by factor of 2
2. Bad initial configuration (overlapping atoms) → Run energy minimization
3. Wrong units (mixing real and LJ) → Check all parameters
4. Force field mismatch → Verify `pair_coeff` values

**Diagnostic:**
```lammps
# Quick test with very small timestep
timestep 0.0001
run 100
```

If this works, gradually increase timestep.

### 10.2 System Doesn't Equilibrate

**Symptom:** Temperature/pressure/density still drifting after 100,000 steps

**Causes:**
1. Equilibration time too short → Run 10× longer
2. Multiple metastable states (glass) → Need enhanced sampling
3. Box too small → Finite-size effects
4. Wrong ensemble (using NVE when you meant NVT)

**Diagnostic:**
Plot running average of property — if slope ≠ 0, not equilibrated.

### 10.3 Nonsensical Results

**Symptom:** Calculated properties don't match literature or expectations

**Checklist:**
- [ ] Units are consistent throughout input file
- [ ] Equilibration data was discarded before averaging
- [ ] Sufficient statistics (>10,000 independent samples)
- [ ] Correct force field for the material
- [ ] Box size large enough (>10× cutoff in each dimension)
- [ ] Temperature/pressure actually controlled (check thermo output)

---

## 11. Validation and Verification

### 11.1 Internal Consistency Checks

**Before trusting any result:**

1. **Energy conservation (NVE test):**
   ```lammps
   fix 1 all nve
   run 100000
   # Check: Total energy drift < 0.01% per ns
   ```

2. **Temperature distribution (NVT test):**
   * Calculate kinetic energy distribution
   * Should follow chi-squared with $N_f$ degrees of freedom
   * Tutorial 7 shows how to check this

3. **Pressure/volume stability (NPT test):**
   * Running average of volume should be flat
   * Pressure should fluctuate around target (±10-20% is normal)

### 11.2 Comparison with Known Results

**For learning/testing:**
* Argon: Well-studied, extensive literature data
* Water (SPC/E, TIP3P): Benchmark properties tabulated
* Proteins: Compare RMSD, Rg with experimental values

**Example: Argon liquid at 94.4 K**

Literature values (from NIST):
* Density: 1.374 g/cm³
* Pressure: 0.1 MPa

Your simulation should match within 1-2%.

---

## 12. From Simulation to Publication

### 12.1 What to Report

**Minimum information for reproducibility:**

1. **System details:**
   * Number of atoms, composition
   * Box size (initial and final if NPT)
   * Initial configuration (lattice, random, from file)

2. **Force field:**
   * Potential type (LJ, AMBER, etc.)
   * All parameters (ε, σ, charges, bonds)
   * Cutoffs and long-range corrections

3. **Simulation protocol:**
   * Timestep
   * Ensemble (NVE/NVT/NPT)
   * Thermostat/barostat type and parameters
   * Equilibration duration
   * Production duration

4. **Analysis method:**
   * How equilibration was detected
   * How averages were calculated
   * Uncertainty estimation method (block averaging)

**Example (Methods section):**
> We performed molecular dynamics simulations of 2000 Argon atoms using the Lennard-Jones potential with ε = 0.238 kcal/mol and σ = 3.405 Å, with a cutoff of 10.0 Å. Simulations were performed using LAMMPS version 20210310 on the Purdue Anvil cluster. The system was equilibrated in the NPT ensemble at T = 94.4 K and P = 1 atm for 100 ps, followed by a 1 ns production run in the NVT ensemble. The Nosé-Hoover thermostat with damping parameter 100 fs was used for temperature control. A timestep of 1 fs was used throughout. Average density was calculated by discarding the first 200 ps of the production run and block-averaging the remaining 800 ps with 10 blocks, yielding ρ = 1.372 ± 0.008 g/cm³ (95% CI).

### 12.2 Module Loading for Reproducibility

**Always document the exact modules used:**

```bash
# Load required modules (Anvil)
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310
```

This ensures others can reproduce your results with the same software environment.

### 12.3 Data Archiving

**Make your simulations reproducible:**
* Save complete input files (GitHub repository)
* Archive final configurations
* Document software version: `lmp -help | head`
* Include analysis scripts (Python notebooks)

---

## 13. Summary: The Simulation Checklist

Before running any production simulation:

**Planning:**
- [ ] Clearly define the research question
- [ ] Choose appropriate system size (balance accuracy vs. cost)
- [ ] Select correct force field for the material
- [ ] Determine required simulation time (literature survey)

**Setup:**
- [ ] Verify all parameters are in consistent units
- [ ] Check initial configuration is reasonable (no overlaps)
- [ ] Test timestep stability (NVE test)
- [ ] Estimate computational cost (wall time × core hours)

**Execution:**
- [ ] Energy minimization (if needed)
- [ ] NVT equilibration (heat to target T)
- [ ] NPT equilibration (if studying density-dependent properties)
- [ ] Production run with sufficient length
- [ ] Monitor job progress (energy, temperature)

**Analysis:**
- [ ] Detect and discard equilibration phase
- [ ] Use block averaging for uncertainties
- [ ] Calculate multiple independent runs for robustness
- [ ] Validate against known results (if available)

**Reporting:**
- [ ] Document all parameters and methods
- [ ] Report uncertainties (mean ± 95% CI)
- [ ] Archive input files and data

---

## 14. Further Reading

### 14.1 Essential References

**Textbooks:**
* **Frenkel & Smit, "Understanding Molecular Simulation"** — Chapter 4 (MD Basics), Chapter 6 (Free Energy)
* **Allen & Tildesley, "Computer Simulation of Liquids"** — The classic reference
* **Tuckerman, "Statistical Mechanics: Theory and Molecular Simulation"** — Rigorous theoretical foundation

**LAMMPS Documentation:**
* **[LAMMPS Manual](https://docs.lammps.org/)** — Comprehensive command reference
* **[LAMMPS Howtos](https://docs.lammps.org/Howto.html)** — Topic-specific guides
* **[LAMMPS Examples](https://github.com/lammps/lammps/tree/develop/examples)** — 100+ example input files

**Papers on Best Practices:**
* **González (2011)** — "Estimation of the relative free energy of solvation from MD" (umbrella sampling)
* **Hess (2002)** — "Determining the shear viscosity of model fluids" (transport properties)
* **Yeh & Hummer (2004)** — "System-size dependence of diffusion coefficients" (finite-size corrections)

### 14.2 Online Resources

* **[LAMMPS on Anvil](https://www.rcac.purdue.edu/knowledge/anvil/software/installing_applications/lammps/provided_module)** — Anvil-specific LAMMPS documentation
* **[Materials Cloud](https://www.materialscloud.org/work/tools/seekpath)** — Visualization and analysis tools
* **[LAMMPS Tutorials (lammpstutorials.github.io)](https://lammpstutorials.github.io/)** — Community tutorials
* **[LAMMPS Workshop Videos](https://www.lammps.org/workshops.html)** — Annual training workshops

:::{note} Class Examples
Sample LAMMPS input files, job scripts, and output data are available in:
```bash
/anvil/projects/x-chm250117/class_examples/
```
These tested examples demonstrate best practices for common simulation tasks.
:::

---

## 15. Looking Ahead

This lecture covered the **core workflow** of running MD simulations. Advanced topics (covered in future lectures or independent study):

* **Enhanced Sampling:** Umbrella sampling, metadynamics, replica exchange
* **Free Energy Calculations:** Thermodynamic integration, Bennett acceptance ratio
* **Advanced Force Fields:** Reactive (ReaxFF), polarizable (Drude), machine-learned (DeepMD)
* **Coarse-Graining:** MARTINI, DPD, dissipative particle dynamics
* **Non-Equilibrium MD:** Shear flows, thermal gradients, shock waves

The skills you've learned — writing input files, running simulations on HPC, analyzing output — form the foundation for all of these advanced techniques.

:::{important} The Bottom Line
Molecular dynamics is a **tool**, not magic. Results are only as good as:
1. Your force field (garbage in = garbage out)
2. Your equilibration (biased sampling = biased results)
3. Your statistics (insufficient sampling = large uncertainties)

Always validate, always check, always question your results.
:::