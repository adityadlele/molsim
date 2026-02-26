# Tutorial 9: Running LAMMPS on Anvil Cluster

**Objective:** Learn to submit LAMMPS jobs to the Anvil supercomputer using SLURM batch scripts.

:::{note} Prerequisites
* Tutorial 2: Accessing Anvil
* Tutorial 3: Linux Survival
* Tutorial 7: LAMMPS Input Files
* Active ACCESS account added to class allocation
:::

---

## 1. Why Use a Job Scheduler?

On your laptop, you click "Run" and the program starts immediately. On a supercomputer shared by hundreds of researchers, this doesn't work. We need a **job scheduler** (also called a **workload manager**) to:

1. **Queue jobs fairly** — First-come, first-served (with priority adjustments)
2. **Allocate resources** — Assign specific nodes and cores to your job
3. **Track usage** — Ensure you don't exceed your allocation quota
4. **Manage conflicts** — Prevent two users from trying to use the same CPU

Anvil uses **SLURM** (Simple Linux Utility for Resource Management), the most widely-used HPC scheduler in the world.

:::{important} The Golden Rule
**Never run heavy computations on the login node.**

When you SSH into Anvil, you land on a "login node" — a shared machine used by everyone to prepare jobs. If you run a LAMMPS simulation there, it will:
* Slow down everyone else's work
* Get automatically killed by the system
* Potentially get your account temporarily suspended

Always use `sbatch` to submit jobs to the compute nodes.
:::

---

## 2. The SLURM Workflow

```text
┌─────────────────┐
│  Your Laptop    │
└────────┬────────┘
         │ SSH
         ↓
┌─────────────────┐
│  Login Node     │  ← You write/edit files here
│  (Shared)       │  ← Never run simulations here!
└────────┬────────┘
         │ sbatch job.sh
         ↓
┌─────────────────┐
│  SLURM Queue    │  ← Job waits for available resources
└────────┬────────┘
         │ Job starts
         ↓
┌─────────────────┐
│  Compute Node   │  ← Simulation runs here
│  (Exclusive)    │  ← You have dedicated cores
└────────┬────────┘
         │ Job completes
         ↓
┌─────────────────┐
│  Output Files   │  ← Results written to $SCRATCH
└─────────────────┘
```

---

## 3. Understanding SLURM Partitions

Anvil divides compute nodes into **partitions** (queues) based on job size and duration.

| Partition | Max Nodes | Max Time | Cores/Node | Use Case |
|-----------|-----------|----------|------------|----------|
| `shared` | 1 (partial) | 48 hours | 128 | Small jobs (testing, short runs) |
| `standard` | 1-25 | 48 hours | 128 | Medium jobs (most production runs) |
| `wide` | 26-100 | 48 hours | 128 | Large parallel jobs |
| `debug` | 1-4 | 30 minutes | 128 | Quick tests (highest priority) |

For this class, you'll primarily use **`shared`** and **`standard`**.

:::{seealso}
Full documentation: [Anvil Queue Policies](https://www.rcac.purdue.edu/knowledge/anvil/run/queues)
:::

---

## 4. Anatomy of a SLURM Batch Script

A SLURM script is a regular bash script with special directives at the top. These directives start with `#SBATCH` and tell SLURM what resources you need.

### 4.1 The Template

```bash
#!/bin/bash
#SBATCH --job-name=argon_md          # Job name (shows in queue)
#SBATCH --account=chm250117          # Class allocation ID
#SBATCH --partition=shared           # Which queue to use
#SBATCH --nodes=1                    # Number of compute nodes
#SBATCH --ntasks-per-node=4          # Number of MPI tasks (cores) per node
#SBATCH --time=02:00:00              # Max runtime (HH:MM:SS)
#SBATCH --output=job_%j.out          # Standard output (%j = job ID)
#SBATCH --error=job_%j.err           # Standard error
#SBATCH --mail-type=END,FAIL         # Email when job ends or fails
#SBATCH --mail-user=your_email@rowan.edu

# ========================================
# SECTION 1: ENVIRONMENT SETUP
# ========================================
module purge                         # Clear any loaded modules
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310

# Print job info to output file
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Working directory: $(pwd)"

# ========================================
# SECTION 2: CHANGE TO WORKING DIRECTORY
# ========================================
cd $SLURM_SUBMIT_DIR

# ========================================
# SECTION 3: RUN THE SIMULATION
# ========================================
lmp -in argon.in

# ========================================
# SECTION 4: CLEANUP AND DIAGNOSTICS
# ========================================
echo "Job completed at: $(date)"
```

### 4.2 Key Directives Explained

| Directive | Purpose | Example |
|-----------|---------|---------|
| `--job-name` | Human-readable name for queue | `argon_500K_nvt` |
| `--account` | Allocation project ID | `chm250117` |
| `--partition` | Which queue to submit to | `shared` or `standard` |
| `--nodes` | Number of compute nodes | `1` (typical for MD) |
| `--ntasks-per-node` | CPU cores to use | `4` to `128` |
| `--time` | Max walltime (auto-killed after) | `02:00:00` = 2 hours |
| `--output` | Where to save stdout | `%j` = job ID |
| `--error` | Where to save stderr | `%j` = job ID |
| `--mail-type` | Email notifications | `BEGIN,END,FAIL` |

:::{warning} Time Limit
If your job exceeds the `--time` limit, it will be **killed immediately** without saving state. Always add a 10-20% buffer. For example, if you estimate 5 hours, request 6 hours.
:::

---

## 5. Step-by-Step: Submitting Your First Job

### 5.1 Prepare Your Files

1. **Log in to Anvil** (via SSH or Open OnDemand terminal)
2. **Navigate to your scratch space:**
   ```bash
   cd $SCRATCH
   mkdir lammps_tutorial
   cd lammps_tutorial
   ```

3. **Copy your LAMMPS input file:**
   ```bash
   nano argon.in
   ```
   Paste the input file from Tutorial 7, then save (`Ctrl+O`, `Enter`, `Ctrl+X`).

### 5.2 Create the Submission Script

```bash
nano submit_argon.sh
```

Paste this script (replace `ACCOUNT_NAME` with your actual allocation):

```bash
#!/bin/bash
#SBATCH --job-name=argon_tutorial
#SBATCH --account=chm250117
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=01:00:00
#SBATCH --output=argon_%j.out
#SBATCH --error=argon_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@rowan.edu

# Load LAMMPS modules
module purge
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310

# Print diagnostics
echo "========================================="
echo "Job started: $(date)"
echo "Running on: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Cores requested: $SLURM_NTASKS"
echo "========================================="

# Run LAMMPS
lmp -in argon.in

echo "========================================="
echo "Job completed: $(date)"
echo "========================================="
```

Save the file.

### 5.3 Submit the Job

```bash
sbatch submit_argon.sh
```

You should see:
```text
Submitted batch job 123456
```

The number `123456` is your **Job ID** — you'll use this to check status.

### 5.4 Monitor the Job

**Check queue status:**
```bash
squeue -u $USER
```

Output:
```text
JOBID    PARTITION  NAME          USER      ST  TIME  NODES  NODELIST
123456   shared     argon_tut     alele     R   0:05  1      a001
```

**Status codes:**
* `PD` = Pending (waiting in queue)
* `R` = Running
* `CG` = Completing (cleaning up)
* No output = Job finished

**Cancel a job:**
```bash
scancel 123456
```

### 5.5 Check the Results

Once the job completes (disappears from `squeue`), check for output files:

```bash
ls -lh
```

You should see:
```text
argon.in
submit_argon.sh
argon_123456.out       # Standard output
argon_123456.err       # Error messages
thermo.out             # LAMMPS thermodynamic data
argon.lammpstrj        # Trajectory file
final_state.data       # Final configuration
```

**View the job log:**
```bash
less argon_123456.out
```

**Check for errors:**
```bash
cat argon_123456.err
```

---

## 6. Parallel Execution with MPI

LAMMPS can use multiple CPU cores via MPI (Message Passing Interface) to run faster.

### 6.1 How Many Cores to Use?

**Rule of thumb:** Use 1 core per ~1000-2000 atoms for optimal performance.

For our 2000-atom system:
* **Good:** 1-4 cores
* **Acceptable:** 8 cores
* **Wasteful:** 32+ cores (communication overhead dominates)

### 6.2 Modified Submission Script for MPI

```bash
#!/bin/bash
#SBATCH --job-name=argon_parallel
#SBATCH --account=chm250117
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16         # Request 16 cores
#SBATCH --time=00:30:00
#SBATCH --output=argon_mpi_%j.out
#SBATCH --error=argon_mpi_%j.err

module purge
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310

echo "Using $SLURM_NTASKS MPI tasks"

# Run LAMMPS with MPI
mpirun -np $SLURM_NTASKS lmp -in argon.in

echo "Simulation completed"
```

**Key change:** `mpirun -np $SLURM_NTASKS lmp -in argon.in`

This launches LAMMPS across multiple cores. LAMMPS automatically divides the simulation box spatially (domain decomposition).

---

## 7. Example Submission Scripts

### 7.1 Production Run (Long Simulation)

```bash
#!/bin/bash
#SBATCH --job-name=argon_production
#SBATCH --account=chm250117
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00              # 24 hours
#SBATCH --output=prod_%j.out
#SBATCH --error=prod_%j.err

module purge
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310

# Important: Run from SCRATCH (fast I/O)
cd $SCRATCH/lammps_production
cp $SLURM_SUBMIT_DIR/argon.in .

# Run simulation
mpirun -np $SLURM_NTASKS lmp -in argon.in

# Copy results back to submission directory
cp thermo.out final_state.data $SLURM_SUBMIT_DIR/

echo "Production run completed at $(date)"
```

### 7.2 Test Run (Quick Debugging)

```bash
#!/bin/bash
#SBATCH --job-name=quick_test
#SBATCH --account=chm250117
#SBATCH --partition=debug            # High priority, short time
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:10:00              # 10 minutes max
#SBATCH --output=test_%j.out

module purge
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310

# Modify input to run only 1000 steps
sed 's/run.*/run 1000/' argon.in > test.in

lmp -in test.in

echo "Test completed"
```

---

## 8. Checking Allocation Usage

To see how many core-hours you've used:

```bash
mybalance
```

Output:
```text
Account: YOUR_ALLOCATION
Service Units Used: 1234
Service Units Available: 100000
```

**Service Units (SUs)** = core-hours. Running on 8 cores for 2 hours = 16 SUs.

---

## 9. Common SLURM Commands

| Command | Purpose | Example |
|---------|---------|---------|
| `sbatch script.sh` | Submit a job | `sbatch run.sh` |
| `squeue -u $USER` | Check your jobs | `squeue -u alele` |
| `scancel JOB_ID` | Cancel a job | `scancel 123456` |
| `scancel -u $USER` | Cancel all your jobs | `scancel -u alele` |
| `sinfo` | Check cluster status | `sinfo -p shared` |
| `scontrol show job JOB_ID` | Detailed job info | `scontrol show job 123456` |
| `sacct -j JOB_ID` | Job accounting info | `sacct -j 123456` |
| `mybalance` | Check allocation usage | `mybalance` |

---

## 10. Troubleshooting Common Errors

### 10.1 "Job Failed to Start"

**Error in `.err` file:**
```text
slurmstepd: error: execve(): lmp: No such file or directory
```

**Fix:** LAMMPS module not loaded. Add `module load lammps` to script.

---

### 10.2 "Permission Denied"

**Error:**
```text
sbatch: error: Batch script contains DOS line breaks
```

**Fix:** File was edited on Windows. Convert to Unix format:
```bash
dos2unix submit_argon.sh
```

---

### 10.3 "Exceeded Time Limit"

**Error in `.out` file:**
```text
slurmstepd: error: *** JOB 123456 ON a001 CANCELLED AT 2026-02-26T10:00:00 DUE TO TIME LIMIT ***
```

**Fix:** Increase `--time` in submission script or reduce simulation steps.

---

### 10.4 "Out of Memory"

**Error:**
```text
slurmstepd: error: Detected 1 oom-kill event(s)
```

**Fix:** 
* Request more memory: add `#SBATCH --mem=16G`
* Reduce system size or trajectory dump frequency

---

## 11. Appendix: Installing Custom LAMMPS Packages

LAMMPS is modular — many features are optional "packages" that must be compiled in. The default Anvil LAMMPS module (version 20210310) includes most common packages, but you may need to install a custom build.

:::{seealso}
Official Anvil LAMMPS Documentation: [LAMMPS on Anvil](https://www.rcac.purdue.edu/knowledge/anvil/software/installing_applications/lammps/provided_module)
:::

### 11.1 Check Installed Packages

```bash
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310
lmp -help | grep -A 50 "Installed packages"
```

This lists all compiled packages (e.g., `MOLECULE`, `KSPACE`, `RIGID`, etc.).

### 11.2 Class Examples Directory

Sample LAMMPS input files and job scripts are available in the shared class directory:

```bash
# List available examples
ls /anvil/projects/x-chm250117/class_examples/

# Copy an example to your workspace
cp /anvil/projects/x-chm250117/class_examples/argon_nvt.in $SCRATCH/
cp /anvil/projects/x-chm250117/class_examples/submit_argon.sh $SCRATCH/
```

These examples are tested and ready to use. You can modify them for your specific needs.

### 11.3 When You Need a Custom Build

You'll need to compile LAMMPS yourself if:
* You need a package not in the default build (e.g., `COLVARS`, `PLUMED`)
* You want GPU acceleration (CUDA/Kokkos)
* You need a development feature from LAMMPS GitHub

### 11.4 Basic Compilation Workflow

**Step 1: Load Build Tools**
```bash
module load gcc/11.2.0 openmpi/4.0.6
module load cmake
```

**Step 2: Download LAMMPS Source**
```bash
cd $SCRATCH
git clone -b stable https://github.com/lammps/lammps.git
cd lammps
```

**Step 3: Create Build Directory**
```bash
mkdir build
cd build
```

**Step 4: Configure with CMake**
```bash
cmake ../cmake \
    -DCMAKE_INSTALL_PREFIX=$HOME/software/lammps \
    -DPKG_MOLECULE=ON \
    -DPKG_KSPACE=ON \
    -DPKG_RIGID=ON \
    -DPKG_COLVARS=ON \
    -DBUILD_MPI=ON
```

**Step 5: Compile**
```bash
make -j 8                    # Use 8 cores
make install                 # Install to $HOME/software/lammps
```

**Step 6: Add to PATH**
Add to your `~/.bashrc`:
```bash
export PATH=$HOME/software/lammps/bin:$PATH
```

Reload:
```bash
source ~/.bashrc
```

**Step 7: Test**
```bash
lmp -help
```

### 11.5 Common Package Options

| Package | Purpose | CMake Flag |
|---------|---------|------------|
| `MOLECULE` | Bonds, angles, molecules | `-DPKG_MOLECULE=ON` |
| `KSPACE` | Long-range electrostatics (PME) | `-DPKG_KSPACE=ON` |
| `RIGID` | Rigid body dynamics | `-DPKG_RIGID=ON` |
| `COLVARS` | Collective variables | `-DPKG_COLVARS=ON` |
| `KOKKOS` | GPU acceleration | `-DPKG_KOKKOS=ON` |
| `USER-REAXFF` | Reactive force fields | `-DPKG_REAXFF=ON` |

### 11.6 Using Your Custom LAMMPS in Jobs

In your submission script, replace:
```bash
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310
```

With:
```bash
module load gcc/11.2.0 openmpi/4.0.6
export PATH=$HOME/software/lammps/bin:$PATH
```

:::{tip} Pre-compiled Containers
For complex builds (especially GPU), consider using Singularity containers. Anvil provides pre-built LAMMPS containers:
```bash
module load singularity
singularity exec /apps/containers/lammps-gpu.sif lmp -in input.in
```
:::

---

## 12. Exercise: Submit a Parameter Sweep

**Task:** Run the same Argon simulation at three different temperatures (300K, 500K, 700K) in parallel.

**Approach:** Use a **job array** to submit 3 jobs at once.

**Script:**
```bash
#!/bin/bash
#SBATCH --job-name=argon_sweep
#SBATCH --account=chm250117
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:30:00
#SBATCH --output=sweep_%A_%a.out          # %A = array job ID, %a = task ID
#SBATCH --array=0-2                       # Run tasks 0, 1, 2

module purge
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310

# Define temperatures for each array index
TEMPS=(2.50 4.17 5.84)                    # T* values for 300K, 500K, 700K
T=${TEMPS[$SLURM_ARRAY_TASK_ID]}

echo "Running temperature $T (Task ID: $SLURM_ARRAY_TASK_ID)"

# Create custom input file for this temperature
sed "s/TEMP_PLACEHOLDER/$T/g" template.in > argon_$T.in

# Run simulation
lmp -in argon_$T.in

echo "Completed T=$T"
```

**Prepare template:**
Create `template.in` where temperature lines use `TEMP_PLACEHOLDER`:
```lammps
velocity        all create TEMP_PLACEHOLDER 87287 dist gaussian
fix             1 all nvt temp TEMP_PLACEHOLDER TEMP_PLACEHOLDER 0.5
```

**Submit:**
```bash
sbatch sweep.sh
```

This launches 3 jobs simultaneously, each with a different temperature.

---

## 13. Summary Checklist

Before submitting a job, verify:

- [ ] Input file is correct and tested
- [ ] Working in `$SCRATCH` (not `$HOME`)
- [ ] `--account` matches your allocation
- [ ] `--time` has 20% buffer
- [ ] Email address is correct
- [ ] Output filenames use `%j` to avoid overwriting
- [ ] LAMMPS module is loaded in script
- [ ] Requested cores match system size (~1 core per 1000-2000 atoms)

---

## 14. Further Reading

* **[Anvil User Guide: Running Jobs](https://www.rcac.purdue.edu/knowledge/anvil/run)** — Official documentation
* **[Anvil LAMMPS Documentation](https://www.rcac.purdue.edu/knowledge/anvil/software/installing_applications/lammps/provided_module)** — LAMMPS-specific guide for Anvil
* **[SLURM Quick Start](https://slurm.schedmd.com/quickstart.html)** — Official SLURM tutorial
* **[Job Arrays Tutorial](https://slurm.schedmd.com/job_array.html)** — Advanced parallel job submission

:::{note} Next Tutorial
**Tutorial 9:** Analyzing LAMMPS output with Python — Learn to process `thermo.out` files and calculate statistical properties.
:::