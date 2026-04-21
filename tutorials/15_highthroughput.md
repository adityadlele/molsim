# Tutorial 15: High-Throughput Simulation Workflows

**Objective:** Automate the generation, submission, and analysis of multiple LAMMPS simulations using Python and SLURM job arrays. Apply these tools to a temperature sweep of copper.

:::{note} Prerequisites
- Tutorials 1–4 (Anvil access, Linux, Python/Conda)
- Lecture 12 (high-throughput concepts)
- Basic LAMMPS familiarity (Tutorials 8–12)
:::

---

## Setup

### Working Directory

All files generated in this tutorial — directories, input scripts, SLURM scripts, and simulation outputs — should live in your **scratch space**, not your home directory. Scratch has far more storage and is the correct location for simulation data.

```bash
cd $SCRATCH
mkdir -p tutorial15
cd tutorial15
```

:::{warning} Jupyter notebooks and working directories
If you open a Jupyter notebook from the Anvil OnDemand launcher, it starts in your **home directory** (`/home/x-yourusername`), not in scratch. This means relative paths like `Path("temperature_sweep")` or `Path("lammps_nvt_template.in")` will resolve to your home directory, not to the tutorial15 folder you just created.

To keep everything in one place, run this at the top of your notebook before any other cells:

```python
import os
os.chdir(f"/anvil/scratch/{os.environ['USER']}/tutorial15")
print(os.getcwd())   # confirm you are in the right place
```

Alternatively, navigate to `$SCRATCH/tutorial15` in the Jupyter file browser and open (or create) your notebook from there.
:::

### Pre-Run Outputs

This tutorial involves submitting real LAMMPS jobs. If jobs are still running during the session, pre-computed outputs are available at:

```
/anvil/projects/x-chm250117/class_examples/tutorial15/temperature_sweep/
```

Copy them to your scratch space to follow the analysis sections immediately:

```bash
cp -r /anvil/projects/x-chm250117/class_examples/tutorial15/temperature_sweep \
      $SCRATCH/tutorial15/prerun
```

---

## The Science: Thermal Properties of Copper

We will run NVT simulations of FCC copper at six temperatures spanning 300 K to 1300 K (just below the experimental melting point of 1358 K). At each temperature we extract the mean potential energy per atom and the mean pressure.

This produces two physically meaningful results:
1. How potential energy increases with temperature (thermal energy storage).
2. How pressure builds at fixed volume as the lattice wants to expand (thermal pressure, related to the thermal expansion coefficient).

---

## Part 0: Bash Basics — Before We Touch the Cluster

The LAMMPS workflows later in this tutorial use the same two bash patterns over and over: loops that build directories, and `sed` substitution that fills in template files. This section introduces both in isolation — no simulation software, no SLURM — so the mechanics are clear before they get buried in a larger script.

---

### Example 1: Generating a Structured Directory Tree

Suppose you are planning a study with three materials, each run at four pressures. You want a clean folder for every combination before a single simulation runs. Doing this by hand for 12 folders is tedious; doing it for 120 is not feasible.

The script below builds the tree automatically and drops a `README` into each folder with the parameters baked in.

```bash
#!/bin/bash
# build_study_tree.sh
# Generates a labeled directory tree for a parameter study.

MATERIALS=(copper aluminum iron)
PRESSURES=(1 10 100 1000)

for mat in "${MATERIALS[@]}"; do
    for P in "${PRESSURES[@]}"; do

        # Build a descriptive folder name
        dir="${mat}/P${P}bar"
        mkdir -p "${dir}"

        # Write a README inside with the run metadata
        cat > "${dir}/README" << EOF
Material : ${mat}
Pressure : ${P} bar
Status   : pending
EOF

        echo "Created: ${dir}"
    done
done

echo "Done. $(ls -d */P* | wc -l) directories created."
```

Save this as `build_study_tree.sh` and run it:

```bash
chmod +x build_study_tree.sh
bash build_study_tree.sh
```

The output tree looks like this:

```
copper/
  P1bar/README
  P10bar/README
  P100bar/README
  P1000bar/README
aluminum/
  P1bar/README
  ...
```

Check one of the README files to confirm the substitution worked:

```bash
cat copper/P100bar/README
```

Expected output:
```
Material : copper
Pressure : 100 bar
Status   : pending
```

A few things to notice. `mkdir -p` creates parent directories as needed — `copper/P1bar/` in one command, no errors if the folder already exists. The `cat > file << EOF ... EOF` pattern (a *here-document*) lets you write a multi-line file inline without a separate template file. The `${mat}` and `${P}` variables are substituted before the content is written, so each README is unique.

---

### Example 2: Template Files and `sed` Substitution

A here-document is convenient for small files, but when the file is long (like a LAMMPS input script with dozens of lines), it is cleaner to keep a **template file** separate and use `sed` to swap in values. This is the pattern used in every workflow later in this tutorial.

Start by creating a template. This one mimics a generic simulation config — the same idea applies directly to LAMMPS input files.

```bash
cat > params_template.conf << 'EOF'
# Simulation configuration
# Generated automatically — do not edit by hand

system      = Cu_FCC
temperature = TEMP_K
pressure    = PRESS_BAR
timestep    = 0.002
run_steps   = 50000
output_dir  = results/T_TEMP_K_P_PRESS_BAR
EOF
```

Now write a loop that generates one filled-in config per temperature and pressure combination:

```bash
#!/bin/bash
# fill_templates.sh

TEMPERATURES=(300 600 900 1200)
PRESSURES=(1 100)

for T in "${TEMPERATURES[@]}"; do
    for P in "${PRESSURES[@]}"; do

        outfile="config_T${T}_P${P}.conf"

        # Chain two sed calls: replace temperature first, then pressure
        sed "s/TEMP_K/${T}/g" params_template.conf \
          | sed "s/PRESS_BAR/${P}/g" \
          > "${outfile}"

        echo "Wrote: ${outfile}"
    done
done
```

Run it:

```bash
chmod +x fill_templates.sh
bash fill_templates.sh
ls *.conf
```

Inspect one output:

```bash
cat config_T600_P100.conf
```

Expected:
```
# Simulation configuration
# Generated automatically — do not edit by hand

system      = Cu_FCC
temperature = 600
pressure    = 100
timestep    = 0.002
run_steps   = 50000
output_dir  = results/T_600_P_100
```

`sed "s/OLD/NEW/g"` replaces every occurrence of `OLD` with `NEW` in the file. Piping two `sed` calls with `|` applies both substitutions in one pass. The template itself is never modified — you can regenerate all configs from scratch at any time by re-running the script.

:::{tip} Combining both patterns
In practice, you will combine Examples 1 and 2: the directory loop creates folders, and `sed` fills in the input file for each one. That is exactly what Part 1 and Part 2 of this tutorial do with real LAMMPS inputs. The only addition is `sbatch` at the end of each loop iteration to hand the job to the cluster.
:::

---

## Part 1: Python-Driven Workflow

### Step 1.1: The LAMMPS Template

Save the following as `$SCRATCH/tutorial15/lammps_nvt_template.in`:

```lammps
# Cu NVT temperature sweep — template
# TEMP_TARGET will be replaced by the Python script

units       metal
atom_style  atomic
boundary    p p p

lattice     fcc 3.615
region      box block 0 5 0 5 0 5
create_box  1 box
create_atoms 1 box

pair_style  eam
pair_coeff  * * /apps/spack/anvil/apps/lammps/20210310-gcc-11.2.0-jzfe7x3/share/lammps/potentials/Cu_u3.eam

# Equilibration: heat to target temperature
velocity    all create TEMP_TARGET 42 dist gaussian
fix         eq all nvt temp TEMP_TARGET TEMP_TARGET 0.1
thermo      500
run         5000
unfix       eq

# Production
reset_timestep 0
fix         prod all nvt temp TEMP_TARGET TEMP_TARGET 0.1
thermo_style custom step temp pe press
thermo      100
run         10000
```

### Step 1.2: Python Workflow Script

Create a Jupyter notebook (or Python script) in your tutorial15 directory:

**Cell 1: Generate input files**

```python
import subprocess
from pathlib import Path

# ── Configuration ─────────────────────────────────────────────────
TEMPERATURES = [300, 500, 700, 900, 1100, 1300]
BASE_DIR     = Path("temperature_sweep")
TEMPLATE     = Path("lammps_nvt_template.in").read_text()

# ── Generate directories and input files ──────────────────────────
for T in TEMPERATURES:
    run_dir = BASE_DIR / f"T{T}"
    run_dir.mkdir(parents=True, exist_ok=True)

    # Replace placeholder with actual temperature
    input_text = TEMPLATE.replace("TEMP_TARGET", str(T))
    (run_dir / "lammps.in").write_text(input_text)

    print(f"Created: {run_dir}/lammps.in")

print(f"\nGenerated {len(TEMPERATURES)} input files.")
```

**Cell 2: Write and submit SLURM scripts**

```python
job_ids = {}

for T in TEMPERATURES:
    run_dir = BASE_DIR / f"T{T}"

    slurm_script = f"""#!/bin/bash
#SBATCH --job-name=cu_T{T}
#SBATCH --account=chm250117
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:20:00
#SBATCH --output=slurm_%j.out

module purge
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310

mpirun -np 4 lmp -in lammps.in > lammps.out
echo "Job complete: T={T} K"
"""
    (run_dir / "submit.sh").write_text(slurm_script)

    result = subprocess.run(
        ["sbatch", "submit.sh"],
        cwd=run_dir,
        capture_output=True,
        text=True
    )

    if result.returncode == 0:
        job_id = result.stdout.strip().split()[-1]
        job_ids[T] = job_id
        print(f"T={T:5d} K  →  job ID {job_id}")
    else:
        print(f"T={T:5d} K  →  SUBMISSION FAILED: {result.stderr.strip()}")

print(f"\nSubmitted {len(job_ids)} jobs.")
```

**Cell 3: Monitor job status**

```python
import subprocess, time

def check_jobs(job_ids):
    """Return a dict of {job_id: status} for our submitted jobs."""
    result = subprocess.run(
        ["squeue", "-u", subprocess.getoutput("whoami"), "-o", "%i %T"],
        capture_output=True, text=True
    )
    running = {}
    for line in result.stdout.strip().split("\n")[1:]:   # skip header
        parts = line.split()
        if len(parts) == 2:
            running[parts[0]] = parts[1]

    status = {}
    for T, jid in job_ids.items():
        status[T] = running.get(jid, "COMPLETED")
    return status

# Check once (run this cell repeatedly to monitor)
status = check_jobs(job_ids)
for T, s in status.items():
    print(f"T={T:5d} K  →  {s}")
```

### Step 1.3: Parse and Analyze Results

Use either your submitted jobs or the pre-run outputs:

```python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

# Point to pre-run outputs if your jobs haven't finished:
# results_dir = Path("prerun/temperature_sweep")
results_dir = Path("temperature_sweep")

def parse_lammps_thermo(filepath, n_skip_fraction=0.2):
    """
    Parse LAMMPS thermo output from a log file.
    Skips the first n_skip_fraction of steps as equilibration.

    Returns dict of {column_name: np.array}
    """
    data = {"Step": [], "Temp": [], "PotEng": [], "Press": []}
    reading = False

    with open(filepath) as f:
        for line in f:
            # Detect thermo header line
            if line.strip().startswith("Step") and "Temp" in line:
                reading = True
                continue
            if reading:
                parts = line.split()
                # Stop at end of run block
                if len(parts) != 4:
                    reading = False
                    continue
                try:
                    data["Step"].append(int(parts[0]))
                    data["Temp"].append(float(parts[1]))
                    data["PotEng"].append(float(parts[2]))
                    data["Press"].append(float(parts[3]))
                except ValueError:
                    reading = False

    # Convert to arrays and skip equilibration
    for k in data:
        data[k] = np.array(data[k])

    n_skip = int(len(data["Step"]) * n_skip_fraction)
    return {k: v[n_skip:] for k, v in data.items()}


# ── Collect results across all temperatures ──────────────────────
records = []

for T in TEMPERATURES:
    log_file = results_dir / f"T{T}" / "lammps.out"
    if not log_file.exists():
        print(f"T={T}: output not found, skipping.")
        continue

    thermo = parse_lammps_thermo(log_file)

    if len(thermo["Temp"]) == 0:
        print(f"T={T}: no thermo data parsed, skipping.")
        continue

    records.append({
        "T_target":   T,
        "T_mean":     np.mean(thermo["Temp"]),
        "T_std":      np.std(thermo["Temp"]),
        "PE_mean":    np.mean(thermo["PotEng"]),
        "PE_std":     np.std(thermo["PotEng"]),
        "Press_mean": np.mean(thermo["Press"]),
        "Press_std":  np.std(thermo["Press"]),
    })
    print(f"T={T:5d} K  →  <T> = {records[-1]['T_mean']:.1f} K  "
          f"<PE> = {records[-1]['PE_mean']:.3f} eV  "
          f"<P> = {records[-1]['Press_mean']:.0f} bar")

df = pd.DataFrame(records)
print(f"\nCollected {len(df)} runs.")
```

**Cell 5: Plot results**

```python
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: potential energy per atom vs temperature
ax = axes[0]
ax.errorbar(df["T_target"], df["PE_mean"], yerr=df["PE_std"],
            fmt='bo-', lw=2, capsize=4, markersize=7)
ax.set_xlabel("Target Temperature (K)", fontsize=12)
ax.set_ylabel("Mean Potential Energy (eV)", fontsize=12)
ax.set_title("Cu: Thermal Energy Storage", fontsize=13, fontweight='bold')
ax.grid(alpha=0.3)

# Right: pressure vs temperature
ax2 = axes[1]
ax2.errorbar(df["T_target"], df["Press_mean"], yerr=df["Press_std"],
             fmt='rs-', lw=2, capsize=4, markersize=7)
ax2.axhline(0, color='gray', ls='--', lw=1)
ax2.set_xlabel("Target Temperature (K)", fontsize=12)
ax2.set_ylabel("Mean Pressure (bar)", fontsize=12)
ax2.set_title("Cu: Thermal Pressure at Fixed Volume", fontsize=13, fontweight='bold')
ax2.grid(alpha=0.3)

plt.suptitle("Copper NVT Temperature Sweep", fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig("cu_temperature_sweep.png", dpi=150)
plt.show()

print("\nINTERPRETATION")
print("="*50)
print("Potential energy increases with T — atoms explore")
print("higher regions of the PES as thermal motion grows.")
print()
print("Pressure increases with T at fixed volume — the")
print("lattice 'wants' to expand but cannot. This thermal")
print("pressure is related to the thermal expansion")
print("coefficient: α = (1/B) × (dP/dT) at const. V.")
```

---

## Part 2: SLURM Job Array

The Python approach above submits one script per run. SLURM job arrays achieve the same result with a single `sbatch` command and are cleaner for large sweeps.

### Step 2.1: The Array Script

Save the following as `$SCRATCH/tutorial15/array_sweep.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=cu_sweep
#SBATCH --account=chm250117
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:20:00
#SBATCH --output=logs/slurm_%A_%a.out
#SBATCH --array=0-5                 # 6 tasks: indices 0 through 5

# ── Parameter array (one value per array index) ───────────────────
TEMPERATURES=(300 500 700 900 1100 1300)
T=${TEMPERATURES[$SLURM_ARRAY_TASK_ID]}

echo "Array task ${SLURM_ARRAY_TASK_ID}: T=${T} K"
echo "Job ID: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

# ── Create run directory ──────────────────────────────────────────
RUN_DIR="array_sweep/T${T}"
mkdir -p "${RUN_DIR}"
cd "${RUN_DIR}"

# ── Generate input file from template ─────────────────────────────
sed "s/TEMP_TARGET/${T}/g" ../../lammps_nvt_template.in > lammps.in

# ── Run simulation ────────────────────────────────────────────────
module purge
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310

mpirun -np 4 lmp -in lammps.in > lammps.out
echo "T=${T} completed at $(date)"
```

### Step 2.2: Submit and Monitor

```bash
# Create log directory
mkdir -p logs

# Make the script executable
chmod +x array_sweep.sh

# Submit the entire sweep in one command
sbatch array_sweep.sh

# Monitor: shows all tasks as cu_sweep_JOBID_0, _1, ..., _5
squeue -u $USER
```

The SLURM output format `%A_%a` writes one log file per task (e.g., `slurm_12345_0.out`), making it easy to identify which task produced which output.

### Step 2.3: Two-Dimensional Sweep

Extending to a 2D parameter space (temperature × system size, for example) requires mapping the 1D array index to a 2D grid:

```bash
#!/bin/bash
#SBATCH --array=0-11        # 6 temperatures × 2 sizes = 12 combinations

TEMPERATURES=(300 500 700 900 1100 1300)
SIZES=(3 5)                 # supercell sizes (3×3×3 and 5×5×5)

N_SIZES=${#SIZES[@]}        # = 2

# Map 1D index to 2D grid
T_IDX=$(( SLURM_ARRAY_TASK_ID / N_SIZES ))
S_IDX=$(( SLURM_ARRAY_TASK_ID % N_SIZES ))

T=${TEMPERATURES[$T_IDX]}
S=${SIZES[$S_IDX]}

RUN_DIR="sweep2d/T${T}_S${S}"
mkdir -p "${RUN_DIR}"
echo "Running T=${T} K, size=${S}×${S}×${S}"
# ... rest of script
```

---

## Part 3: Storing Results with ASE Database

For ongoing projects, storing results in a queryable database is cleaner than CSV files. ASE provides a lightweight database that requires no server setup.

```python
from ase.db import connect
from pathlib import Path
import numpy as np

# Create (or open) a database file
db = connect("cu_sweep_results.db")

# Store results from our parsed DataFrame
for _, row in df.iterrows():
    # Write metadata — no atoms object needed for pure result storage
    # Use a dummy atoms object or store as key-value pairs
    db.write(
        None,                          # no atoms structure needed
        T_target   = int(row.T_target),
        T_mean     = float(row.T_mean),
        T_std      = float(row.T_std),
        PE_mean    = float(row.PE_mean),
        PE_std     = float(row.PE_std),
        Press_mean = float(row.Press_mean),
        potential  = "EAM_Cu_u3",
        system     = "Cu_FCC_5x5x5",
    )

print(f"Stored {len(df)} entries in cu_sweep_results.db")

# ── Query the database ─────────────────────────────────────────────
print("\nAll runs above 800 K:")
for row in db.select("T_target>800"):
    print(f"  T={row.T_target} K  <PE>={row.PE_mean:.4f} eV  "
          f"<P>={row.Press_mean:.0f} bar")

# Query and plot
import matplotlib.pyplot as plt

T_vals  = []
PE_vals = []
for row in db.select(potential="EAM_Cu_u3"):
    T_vals.append(row.T_target)
    PE_vals.append(row.PE_mean)

# Sort by temperature
idx    = np.argsort(T_vals)
T_vals = np.array(T_vals)[idx]
PE_vals = np.array(PE_vals)[idx]

plt.figure(figsize=(8, 5))
plt.plot(T_vals, PE_vals, 'bo-', lw=2, markersize=8)
plt.xlabel("Temperature (K)", fontsize=12)
plt.ylabel("Mean PE (eV)", fontsize=12)
plt.title("Query from ASE database", fontsize=12)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("from_database.png", dpi=150)
plt.show()
```

The database file `cu_sweep_results.db` is a single portable file. Query it from the command line with:

```bash
ase db cu_sweep_results.db
ase db cu_sweep_results.db T_target=900
```

---

## Summary

| Approach | Best for | Submit command |
|:---------|:---------|:---------------|
| Python loop + subprocess | Full control, adaptive logic | One `sbatch` per run |
| SLURM job array | Clean large sweeps, same parameters | Single `sbatch` |
| Snakemake | Multi-step pipelines with dependencies | `snakemake --executor slurm` |

**Key habits for high-throughput work:**
- Always use a template input file — never edit copies by hand.
- Log both the job ID and the parameters together so you can trace results back to inputs.
- Check a few output files manually before trusting your parser on all 100.
- Store results in a structured format (CSV, ASE database, HDF5) from the beginning — not as raw log files.
