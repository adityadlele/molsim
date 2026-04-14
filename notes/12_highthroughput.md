# High-Throughput Molecular Simulations

---

## 1. Why High-Throughput?

A single MD simulation answers one question: "What happens to this specific system under these specific conditions?" But real research and engineering problems rarely have just one question.

Consider these common scenarios:

- **Parameter screening:** How does the diffusion coefficient of argon change across temperatures from 100 K to 1000 K?
- **Materials discovery:** Which of 50 candidate alloy compositions has the lowest thermal resistance?
- **Convergence testing:** How many atoms do I need before bulk properties stop changing with system size?
- **Statistical sampling:** My result has large variance — I need 10 independent runs with different random seeds to get a reliable mean.

Running these one at a time — manually editing input files, waiting for each job, copying outputs — is untenable. A study with 50 conditions × 5 replicates = 250 simulations would take weeks of human time even if each simulation runs in an hour.

**High-throughput molecular simulation** is the practice of automating the generation, submission, monitoring, and analysis of large numbers of simulations. The goal is to make the computer do the repetitive work so you can focus on interpreting results.

---

## 2. The Anatomy of a Simulation Workflow

Before automating, it helps to identify the discrete steps in any simulation workflow:

```
1. PARAMETERIZATION
   └── Generate input files (structure, force field, conditions)

2. SUBMISSION
   └── Submit job to scheduler (SLURM on Anvil)

3. EXECUTION
   └── Simulation runs on compute nodes

4. MONITORING
   └── Check job status, detect failures

5. POST-PROCESSING
   └── Parse output files, extract properties

6. AGGREGATION
   └── Collect results across all runs, build dataset

7. ANALYSIS
   └── Plot, fit, interpret
```

Each step can be automated. Steps 1, 5, and 6 are pure file and string manipulation — Python handles these naturally. Steps 2–4 interface with SLURM. Step 7 is your science.

The key insight: once you write a workflow that handles one simulation correctly, running 100 is trivial.

---

## 3. Bash Scripting for Simple Workflows

For straightforward parameter sweeps, a bash script is often the fastest tool to reach for. Bash has native support for loops, variable substitution, and calling SLURM commands.

### 3.1 Template-Based Input Generation

The standard approach is to maintain a **template input file** with placeholder variables, then use `sed` (stream editor) to substitute values:

**Template: `lammps_template.in`**
```bash
units        metal
atom_style   atomic
boundary     p p p

lattice      fcc 3.615
region       box block 0 4 0 4 0 4
create_box   1 box
create_atoms 1 box

pair_style   eam
pair_coeff   * * Cu_u3.eam

velocity     all create TEMPERATURE_PLACEHOLDER 12345

fix          1 all nvt temp TEMPERATURE_PLACEHOLDER TEMPERATURE_PLACEHOLDER 0.1
thermo       100
run          10000
```

**Submission loop: `submit_sweep.sh`**
```bash
#!/bin/bash

TEMPERATURES=(300 500 700 900 1100)

for T in "${TEMPERATURES[@]}"; do

    # Create a directory for this run
    mkdir -p run_T${T}
    cd run_T${T}

    # Substitute temperature into the template
    sed "s/TEMPERATURE_PLACEHOLDER/${T}/g" ../lammps_template.in > lammps.in

    # Write a SLURM submit script for this run
    cat > submit.sh << EOF
#!/bin/bash
#SBATCH --job-name=cu_T${T}
#SBATCH --account=chm250117
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:30:00
#SBATCH --output=slurm_%j.out

module use /anvil/projects/x-chm250117/etc/modules
module load lammps

mpirun -np 4 lmp -in lammps.in > lammps.out
EOF

    sbatch submit.sh
    cd ..

done
echo "Submitted ${#TEMPERATURES[@]} jobs."
```

Running `bash submit_sweep.sh` generates five directories, each with its own input file and SLURM script, and submits them all. The total human interaction time is under one second.

### 3.2 Monitoring Jobs

```bash
# Check status of all running jobs
squeue -u $USER

# Watch jobs update every 5 seconds
watch -n 5 squeue -u $USER

# Cancel all your jobs if something goes wrong
scancel -u $USER
```

### 3.3 Limitations of Bash

Bash works well for simple loops but becomes unwieldy when:
- Logic depends on simulation output (adaptive workflows)
- You need to handle failures and resubmissions gracefully
- You want to track provenance (which input produced which output)
- The parameter space is more than one dimension

For these cases, Python is the better tool.

---

## 4. Python-Driven Workflows

Python gives you the full power of a programming language for workflow automation: data structures, error handling, logging, file parsing, and scientific libraries all in one environment.

### 4.1 Core Tools

**`pathlib`:** Modern path manipulation (replacing `os.path`).

**`subprocess`:** Run shell commands from Python, including `sbatch`.

**`string.Template` or f-strings:** Input file generation.

**`pandas`:** Collecting results into a structured table.

### 4.2 A Python Workflow Example

```python
import subprocess
from pathlib import Path
import time

# ── Parameters ────────────────────────────────────────────────────
temperatures = [300, 500, 700, 900, 1100]
base_dir     = Path("cu_sweep")
template     = Path("lammps_template.in").read_text()

# ── Generate and submit ───────────────────────────────────────────
job_ids = {}

for T in temperatures:
    run_dir = base_dir / f"T{T}"
    run_dir.mkdir(parents=True, exist_ok=True)

    # Write input file
    input_text = template.replace("TEMPERATURE_PLACEHOLDER", str(T))
    (run_dir / "lammps.in").write_text(input_text)

    # Write SLURM script
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name=cu_T{T}
#SBATCH --account=chm250117
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:30:00
#SBATCH --output=slurm_%j.out

module use /anvil/projects/x-chm250117/etc/modules
module load lammps
mpirun -np 4 lmp -in lammps.in > lammps.out
"""
    (run_dir / "submit.sh").write_text(slurm_script)

    # Submit and capture job ID
    result = subprocess.run(
        ["sbatch", "submit.sh"],
        cwd=run_dir,
        capture_output=True,
        text=True
    )
    job_id = result.stdout.strip().split()[-1]
    job_ids[T] = job_id
    print(f"T={T} K  →  job {job_id}")

print(f"\nSubmitted {len(job_ids)} jobs.")
```

### 4.3 Parsing Results

Once jobs finish, collect results with a parsing loop:

```python
import pandas as pd
import numpy as np

results = []

for T in temperatures:
    log_file = base_dir / f"T{T}" / "lammps.out"
    
    if not log_file.exists():
        print(f"T={T}: output not found, skipping.")
        continue
    
    # Parse LAMMPS thermo output
    temps_inst, pressures = [], []
    reading = False
    
    with open(log_file) as f:
        for line in f:
            if line.startswith("Step"):
                reading = True
                continue
            if reading:
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        temps_inst.append(float(parts[1]))
                        pressures.append(float(parts[2]))
                    except ValueError:
                        reading = False
    
    if temps_inst:
        # Take mean over production steps (skip first 20% as equilibration)
        n_equil = len(temps_inst) // 5
        results.append({
            "T_target": T,
            "T_mean":   np.mean(temps_inst[n_equil:]),
            "P_mean":   np.mean(pressures[n_equil:]),
        })

df = pd.DataFrame(results)
print(df)
df.to_csv("sweep_results.csv", index=False)
```

---

## 5. SLURM Job Arrays

For large parameter sweeps, SLURM **job arrays** are cleaner than submitting individual scripts. A single `sbatch` command launches many jobs, each identified by an array index `$SLURM_ARRAY_TASK_ID`.

### 5.1 Basic Job Array Script

```bash
#!/bin/bash
#SBATCH --job-name=cu_sweep
#SBATCH --account=chm250117
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:30:00
#SBATCH --output=logs/slurm_%A_%a.out   # %A = job ID, %a = array index
#SBATCH --array=0-4                     # 5 tasks: indices 0,1,2,3,4

# Define parameter values (one per array index)
TEMPERATURES=(300 500 700 900 1100)
T=${TEMPERATURES[$SLURM_ARRAY_TASK_ID]}

# Create and enter run directory
mkdir -p run_T${T}
cd run_T${T}

# Generate input file from template
sed "s/TEMPERATURE_PLACEHOLDER/${T}/g" ../lammps_template.in > lammps.in

# Load modules and run
module use /anvil/projects/x-chm250117/etc/modules
module load lammps
mpirun -np 4 lmp -in lammps.in > lammps.out

echo "T=${T} completed."
```

Submit with:
```bash
mkdir -p logs
sbatch array_sweep.sh
```

SLURM creates 5 jobs with IDs like `12345_0`, `12345_1`, ..., `12345_4`. They all run concurrently (subject to queue availability).

### 5.2 Controlling Concurrency

```bash
#SBATCH --array=0-99%10    # 100 jobs, at most 10 running at once
```

This prevents flooding the queue and is good practice on shared resources.

### 5.3 Two-Dimensional Parameter Sweeps

For a grid of parameters (e.g., temperature × pressure), map the 1D array index to 2D coordinates:

```bash
TEMPERATURES=(300 500 700)
PRESSURES=(1 10 100)
N_P=${#PRESSURES[@]}   # = 3

T_IDX=$(( SLURM_ARRAY_TASK_ID / N_P ))
P_IDX=$(( SLURM_ARRAY_TASK_ID % N_P ))

T=${TEMPERATURES[$T_IDX]}
P=${PRESSURES[$P_IDX]}
```

A `--array=0-8` covers all 9 combinations.

---

## 6. Workflow Managers

For research-grade workflows with dependencies between steps, dedicated workflow managers provide job tracking, failure recovery, and provenance recording that bash/Python loops cannot easily match.

### 6.1 Snakemake

[Snakemake](https://snakemake.readthedocs.io) is a Python-based workflow manager that defines rules connecting input files to output files. It automatically determines which jobs need to run (or rerun after failure) based on file timestamps.

```python
# Snakefile
TEMPERATURES = [300, 500, 700, 900, 1100]

rule all:
    input:
        expand("results/T{T}/msd.dat", T=TEMPERATURES)

rule run_md:
    input:  "templates/lammps.in"
    output: "results/T{T}/lammps.out"
    shell:
        """
        mkdir -p results/T{wildcards.T}
        sed 's/TEMP/{wildcards.T}/g' {input} > results/T{wildcards.T}/lammps.in
        lmp -in results/T{wildcards.T}/lammps.in > {output}
        """

rule compute_msd:
    input:  "results/T{T}/lammps.out"
    output: "results/T{T}/msd.dat"
    script: "scripts/compute_msd.py"
```

Run with `snakemake --cores 4`. Snakemake handles dependency resolution automatically and can submit jobs to SLURM with `--executor slurm`.

### 6.2 AiiDA

[AiiDA](https://www.aiida.net) is a more heavyweight framework designed specifically for computational science. It stores every calculation, input, and output in a database (full provenance). It is widely used in the DFT community and increasingly in ML potential workflows.

**Use AiiDA when:**
- You need a permanent record of every calculation (publication-grade provenance)
- Workflows span thousands of calculations across multiple HPC systems
- You are contributing to a shared database (e.g., Materials Cloud)

**Use Snakemake or plain Python when:**
- Your workflow is on one cluster
- You need something running in an afternoon
- Students or collaborators need to read and modify the workflow

### 6.3 ASE Database

For smaller-scale studies, ASE's built-in database provides lightweight result storage without a full workflow manager:

```python
from ase.db import connect

db = connect('simulations.db')

# Store a result
db.write(atoms, temperature=300, diffusion_coeff=2.4e-9, converged=True)

# Query results
for row in db.select('temperature>500,converged=True'):
    print(row.temperature, row.diffusion_coeff)
```

The database is a single file, portable, and queryable from Python or the command line (`ase db simulations.db`).

---

## 7. GPU-Accelerated MD

### 7.1 Why GPUs?

A modern CPU core executes instructions sequentially and has perhaps 8–32 cores per chip. A GPU (Graphics Processing Unit) has thousands of simpler cores designed for massively parallel arithmetic — exactly the kind of computation needed for force evaluation across millions of atom pairs.

For MD simulations, GPU acceleration typically delivers:
- **5–50× speedup** over CPU-only LAMMPS for typical biomolecular or materials systems
- Near-linear scaling with GPU count for large systems
- Access to simulation timescales (microseconds) and system sizes (tens of millions of atoms) not practical on CPUs

### 7.2 GPU Support in LAMMPS

LAMMPS has three GPU-related packages:

**GPU package:** Offloads pair force calculations to the GPU. Works with CUDA (NVIDIA) and OpenCL. Best for medium systems (~10k–1M atoms). Enable with:
```bash
lmp -sf gpu -pk gpu 1 -in lammps.in
```

**KOKKOS package:** A performance-portable backend supporting CUDA, HIP (AMD), and OpenMP. Better for very large systems and modern GPU architectures:
```bash
lmp -k on g 1 -sf kk -pk kokkos cuda/aware yes -in lammps.in
```

**Intel package:** Optimizes for Intel CPUs and Xeon Phi accelerators (less relevant now).

### 7.3 GPU Jobs on Anvil

Anvil has dedicated GPU nodes. A typical GPU SLURM script:

```bash
#!/bin/bash
#SBATCH --job-name=lammps_gpu
#SBATCH --account=chm250117
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --time=02:00:00
#SBATCH --output=slurm_%j.out

module use /anvil/projects/x-chm250117/etc/modules
module load lammps/gpu

lmp -sf gpu -pk gpu 1 -in lammps.in > lammps.out
```

### 7.4 When Does GPU Help?

GPU acceleration is not always faster. The GPU must be kept busy enough to offset the overhead of transferring data between CPU and GPU memory.

| System Size | Recommendation |
|:------------|:---------------|
| < 5,000 atoms | CPU often faster (GPU overhead dominates) |
| 5,000–100,000 atoms | GPU typically 5–20× faster |
| > 100,000 atoms | GPU strongly recommended; multiple GPUs scale well |

**Pair style matters too.** Lennard-Jones and EAM have excellent GPU implementations. ReaxFF GPU support exists but is less mature. Classical force fields (AMBER, CHARMM) on GPU are well-supported and this is the standard for biomolecular production simulations.

### 7.5 ML Potentials on GPU

MACE and other ML potentials run significantly faster on GPU. In the course environment, MACE runs on CPU, but in research settings GPU inference can deliver 10–50× speedups:

```python
# MACE on GPU (if available)
from mace.calculators import mace_mp
calc = mace_mp(model="medium", device="cuda", default_dtype="float32")
```

For production ML potential MD simulations, using a GPU is essentially mandatory for any system beyond a few hundred atoms.

### 7.6 GROMACS: A GPU-First Code

[GROMACS](https://www.gromacs.org) is worth mentioning here because it was designed with GPU acceleration as a first-class concern. For biomolecular simulations (proteins, lipids, nucleic acids with explicit solvent), GROMACS on GPU is often 100× faster than LAMMPS on CPU for the same system. If your research involves biomolecules rather than materials, GROMACS is the standard tool.

---

## 8. Summary: Choosing Your Approach

| Scale | Approach | Tool |
|:------|:---------|:-----|
| 2–10 runs | Manual or simple loop | Bash `for` loop |
| 10–100 runs | Parameter sweep | Python + `subprocess`, or SLURM job array |
| 100–10,000 runs | Automated workflow | Snakemake + SLURM |
| 10,000+ runs | Full workflow manager | AiiDA, FireWorks |
| Any scale, large systems | Hardware acceleration | LAMMPS GPU/KOKKOS package |
| Biomolecules at scale | GPU-native code | GROMACS |

**Practical advice for this course:** For assignments and the final project, Python-driven submission with SLURM job arrays covers everything you will need. Workflow managers are worth learning if you continue into research — AiiDA in particular is becoming standard in the computational materials community.

---

## Further Reading

- LAMMPS GPU package documentation: [https://docs.lammps.org/Speed_gpu.html](https://docs.lammps.org/Speed_gpu.html)
- Snakemake documentation: [https://snakemake.readthedocs.io](https://snakemake.readthedocs.io)
- AiiDA documentation: [https://aiida.readthedocs.io](https://aiida.readthedocs.io)
- Anvil GPU resources: [https://www.rcac.purdue.edu/compute/anvil](https://www.rcac.purdue.edu/compute/anvil)
