# Tutorial 9: Running LAMMPS on Anvil Cluster

**Objective:** Learn to submit LAMMPS jobs to the Anvil supercomputer using SLURM batch scripts.

:::{note} Prerequisites
* Tutorial 2: Accessing Anvil
* Tutorial 7: LAMMPS Input Files (you should have `argon.in` ready)
* Active ACCESS account added to class allocation
:::

:::{important} In-Class Lab Activity
This is a hands-on lab. We'll submit an actual job to the cluster together. Follow each step carefully and ask for help if stuck!
:::

---

## 1. The Problem: Why Can't I Just Run LAMMPS?

You might think: "I have `argon.in`, I loaded LAMMPS, why can't I just type `lmp -in argon.in` and walk away?"

**Answer:** The login node is shared by everyone. If you run heavy calculations there:
* ❌ You'll slow down everyone else
* ❌ Your job will get killed automatically
* ❌ Your account might get suspended

**Solution:** Use the **job scheduler (SLURM)** to run on dedicated compute nodes.

---

## 2. How the Cluster Works

```text
┌─────────────────┐
│  Your Laptop    │
└────────┬────────┘
         │ SSH
         ↓
┌─────────────────┐
│  Login Node     │  ← You are here (writing files)
│  (Shared)       │  ← DON'T run simulations here!
└────────┬────────┘
         │ sbatch script.sh
         ↓
┌─────────────────┐
│  SLURM Queue    │  ← Your job waits its turn
└────────┬────────┘
         │ Resources available
         ↓
┌─────────────────┐
│  Compute Node   │  ← Simulation runs here
│  (Your cores)   │  ← You have exclusive access
└────────┬────────┘
         │ Job completes
         ↓
┌─────────────────┐
│  Output Files   │  ← Results appear in your directory
└─────────────────┘
```

---

## 3. Step-by-Step: Submitting Jobs from Your Folder Structure

In Tutorial 7, you created organized folders with 5 different simulations. Now we'll submit them as jobs.

### Overview: What We're Submitting

```text
lammps_tutorial/
├── NVT_simulations/
│   ├── T_300K/     → Job 1
│   ├── T_500K/     → Job 2
│   └── T_700K/     → Job 3
└── NPT_simulations/
    ├── P_1atm/     → Job 4
    └── P_10atm/    → Job 5
```

We'll submit **5 separate jobs** (one per folder). Each runs independently on different compute nodes.

---

## 4. Creating Submission Scripts

### Strategy: One Script Per Simulation

Each folder needs its own submission script that:
1. Navigates to the correct folder
2. Runs the specific input file
3. Saves output in that folder

### Step 1: Submit the First NVT Job (T=300K)

```bash
cd $SCRATCH/lammps_tutorial/NVT_simulations/T_300K
nano submit.sh
```

**Paste this script:**

```bash
#!/bin/bash
#SBATCH --job-name=NVT_300K
#SBATCH --account=chm250117
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:30:00
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

echo "========================================="
echo "Job: NVT at 300K"
echo "Started: $(date)"
echo "Directory: $(pwd)"
echo "========================================="

# Load LAMMPS
module purge
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310

# Run simulation
lmp -in argon_300K.in

echo "========================================="
echo "Completed: $(date)"
echo "========================================="
```

Save: `Ctrl+O`, `Enter`, `Ctrl+X`

**Submit the job:**
```bash
sbatch submit.sh
```

You'll see:
```text
Submitted batch job 123456
```

**Check it's running:**
```bash
squeue -u $USER
```

### Step 2: Submit Remaining NVT Jobs

**For T=500K:**
```bash
cd ../T_500K
nano submit.sh
```

Same script, but change:
```bash
#SBATCH --job-name=NVT_500K
echo "Job: NVT at 500K"
lmp -in argon_500K.in
```

Submit:
```bash
sbatch submit.sh
```

**For T=700K:**
```bash
cd ../T_700K
nano submit.sh
```

Change:
```bash
#SBATCH --job-name=NVT_700K
echo "Job: NVT at 700K"
lmp -in argon_700K.in
```

Submit:
```bash
sbatch submit.sh
```

### Step 3: Submit NPT Jobs

**For P=1atm:**
```bash
cd ../../NPT_simulations/P_1atm
nano submit.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=NPT_1atm
#SBATCH --account=chm250117
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=01:00:00
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

echo "========================================="
echo "Job: NPT at 1 atm"
echo "Started: $(date)"
echo "========================================="

module purge
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310

lmp -in argon_npt_1atm.in

echo "========================================="
echo "Completed: $(date)"
echo "========================================="
```

Submit:
```bash
sbatch submit.sh
```

**For P=10atm:**
```bash
cd ../P_10atm
nano submit.sh
```

Change:
```bash
#SBATCH --job-name=NPT_10atm
echo "Job: NPT at 10 atm"
lmp -in argon_npt_10atm.in
```

Submit:
```bash
sbatch submit.sh
```

---

## 5. Monitoring All Your Jobs

**Check all 5 jobs:**
```bash
squeue -u $USER
```

You should see:
```text
JOBID    PARTITION  NAME       USER      ST  TIME  NODES
123456   shared     NVT_300K   username  R   0:02  1
123457   shared     NVT_500K   username  R   0:01  1
123458   shared     NVT_700K   username  PD  0:00  1
123459   shared     NPT_1atm   username  R   0:01  1
123460   shared     NPT_10atm  username  PD  0:00  1
```

**Status meanings:**
* `R` = Running
* `PD` = Pending (waiting for resources)

**Watch live updates:**
```bash
watch -n 5 squeue -u $USER
```

Press `Ctrl+C` to stop.

:::{tip} Cluster Limits
The `shared` partition allows multiple small jobs simultaneously. All 5 might run at once, or some might wait if the cluster is busy.
:::

---

## 6. Collecting Results

After all jobs complete (5-15 minutes), check each folder for output.

### Step 1: Quick Check All Folders

```bash
cd $SCRATCH/lammps_tutorial

# Check NVT results
ls NVT_simulations/T_300K/
ls NVT_simulations/T_500K/
ls NVT_simulations/T_700K/

# Check NPT results
ls NPT_simulations/P_1atm/
ls NPT_simulations/P_10atm/
```

Each folder should now contain:
```text
submit.sh
argon_*.in
job_*.out
job_*.err
thermo_*.out         ← Temperature/energy data
argon_*.lammpstrj    ← Trajectory
final_*.data         ← Final state
log.lammps
```

### Step 2: Verify Successful Completion

**Check job outputs:**
```bash
cd NVT_simulations/T_300K
tail -5 job_*.out
```

Look for:
```text
=========================================
Completed: [timestamp]
=========================================
```

**Check for errors:**
```bash
cat job_*.err
```

Should be empty. If not, read the error message!

### Step 3: Quick Look at Results

**Temperature data for 300K:**
```bash
head -20 thermo_300K.out
```

You should see temperatures fluctuating around 2.50 (target).

**Compare all three NVT temperatures:**
```bash
cd $SCRATCH/lammps_tutorial/NVT_simulations

# Extract average temperatures (after equilibration)
tail -100 T_300K/thermo_300K.out | awk '{sum+=$2; count++} END {print "T_300K avg:", sum/count}'
tail -100 T_500K/thermo_500K.out | awk '{sum+=$2; count++} END {print "T_500K avg:", sum/count}'
tail -100 T_700K/thermo_700K.out | awk '{sum+=$2; count++} END {print "T_700K avg:", sum/count}'
```

:::{note} What You Should See
The average temperatures should be close to your targets:
* T_300K: ~2.50
* T_500K: ~4.17
* T_700K: ~5.84

Small deviations (±0.05) are normal due to thermal fluctuations.
:::

---

## 7. Advanced: Using a Loop to Submit Multiple Jobs

**Pro tip:** Instead of submitting each job manually, use a bash loop.

### Create All NVT Submission Scripts at Once

```bash
cd $SCRATCH/lammps_tutorial/NVT_simulations

# Loop through temperature folders
for temp_dir in T_300K T_500K T_700K; do
    cd $temp_dir
    
    # Extract temperature from folder name (e.g., T_300K → 300K)
    temp_label=${temp_dir#T_}
    
    # Create submission script
    cat > submit.sh << EOF
#!/bin/bash
#SBATCH --job-name=NVT_${temp_label}
#SBATCH --account=chm250117
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:30:00
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

echo "Job: NVT at ${temp_label}"
echo "Started: \$(date)"

module purge
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310

lmp -in argon_${temp_label}.in

echo "Completed: \$(date)"
EOF
    
    # Submit the job
    sbatch submit.sh
    
    # Go back to parent directory
    cd ..
done
```

This creates and submits all 3 NVT jobs with one command!

---

## 8. What If a Job Failed?

### Diagnosing Failures

**Symptom 1: Job completes in < 5 seconds**

Check error file:
```bash
cat job_*.err
```

**Common error:** "Cannot open input script"
```text
Fix: Check filename matches exactly
ls -lh argon_*.in
Make sure the filename in submit.sh matches
```

**Common error:** "lmp: command not found"
```text
Fix: Module loading failed
Add these lines to your submit.sh:
module purge
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310
```

**Symptom 2: Job runs but no thermo output**

Check LAMMPS log:
```bash
tail -50 log.lammps
```

Look for ERROR messages.

**Symptom 3: Job killed before completion**

```bash
tail job_*.out
```

Look for:
```text
slurmstepd: error: *** JOB CANCELLED AT [time] DUE TO TIME LIMIT ***
```

**Fix:** Increase time in `#SBATCH --time=01:00:00`

### Resubmitting Failed Jobs

```bash
# Fix the problem (edit submit.sh or input file)
nano submit.sh

# Resubmit
sbatch submit.sh
```

Old job output files won't be overwritten (they have unique job IDs).

---

## 9. Useful SLURM Commands

| Command | Purpose |
|---------|---------|
| `squeue -u $USER` | Check your jobs |
| `squeue -u $USER --start` | Estimated start times |
| `scancel 123456` | Cancel specific job |
| `scancel -u $USER` | Cancel ALL your jobs |
| `scontrol show job 123456` | Detailed job info |
| `sacct -j 123456` | Job history (after completion) |

---

## 10. Summary Checklist

After this tutorial, you should have:

- [ ] 5 folders each with a `submit.sh` script
- [ ] Successfully submitted 5 jobs (3 NVT + 2 NPT)
- [ ] All jobs completed without errors
- [ ] `thermo_*.out` files in each folder
- [ ] Understand how to check job status with `squeue`
- [ ] Know how to diagnose job failures

**File structure check:**
```bash
cd $SCRATCH/lammps_tutorial
find . -name "thermo_*.out" | wc -l
```

Should output: `5` (one per simulation)

**Next Tutorial:** Analyzing these results with Python (Tutorial 9)

---

## 11. Quick Reference

**Submit job from any folder:**
```bash
sbatch submit.sh
```

**Check all your running jobs:**
```bash
squeue -u $USER
```

**Cancel a stuck job:**
```bash
scancel JOB_ID
```

**Check if job finished successfully:**
```bash
tail -5 job_*.out
# Should see "Completed: [timestamp]"
```

**View real-time progress:**
```bash
tail -f thermo_*.out
# Press Ctrl+C to stop
```

---

## 12. Class Examples

Complete working examples available:
```bash
ls /anvil/projects/x-chm250117/class_examples/
```

If you need to start fresh:
```bash
cd $SCRATCH
rm -rf lammps_tutorial
cp -r /anvil/projects/x-chm250117/class_examples/lammps_tutorial ./
```

---

## 13. Further Reading

* **[Anvil User Guide](https://www.rcac.purdue.edu/knowledge/anvil/run)** — Complete SLURM documentation
* **[Anvil LAMMPS](https://www.rcac.purdue.edu/knowledge/anvil/software/installing_applications/lammps/provided_module)** — LAMMPS on Anvil
* **[SLURM Commands](https://slurm.schedmd.com/pdfs/summary.pdf)** — Quick reference PDF

Let's submit the Argon simulation from Tutorial 7.

### Step 1: Make Sure You Have Your Input File

```bash
cd $SCRATCH/lammps_tutorial
ls -lh argon.in
```

You should see your `argon.in` file. If not, go back to Tutorial 7!

### Step 2: Create a Submission Script

A **submission script** is a text file that tells SLURM:
* What resources you need (cores, time, memory)
* What commands to run

Create the script:
```bash
nano submit_argon.sh
```

### Step 3: Copy This Script

**Paste the following into nano** (carefully!):

```bash
#!/bin/bash
#SBATCH --job-name=my_first_job
#SBATCH --account=chm250117
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:30:00
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

# Print job info
echo "========================================="
echo "Job started: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "========================================="

# Load LAMMPS
module purge
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310

# Run the simulation
lmp -in argon.in

echo "========================================="
echo "Job finished: $(date)"
echo "========================================="
```

**Save and exit:** `Ctrl+O`, `Enter`, `Ctrl+X`

### Step 4: Understand What You Just Wrote

The lines starting with `#SBATCH` are **directives** (instructions to SLURM):

| Line | Meaning |
|------|---------|
| `--job-name=my_first_job` | Name shown in queue |
| `--account=chm250117` | Class allocation ID |
| `--partition=shared` | Which queue (shared = small jobs) |
| `--nodes=1` | Use 1 compute node |
| `--ntasks-per-node=4` | Use 4 CPU cores |
| `--time=00:30:00` | Max 30 minutes (HH:MM:SS) |
| `--output=job_%j.out` | Save output (

%j = job ID) |
| `--error=job_%j.err` | Save errors separately |

### Step 5: Submit Your Job

```bash
sbatch submit_argon.sh
```

You'll see:
```text
Submitted batch job 123456
```

**That number is your Job ID!** Write it down.

### Step 6: Check If It's Running

```bash
squeue -u $USER
```

Output:
```text
JOBID    PARTITION  NAME          USER      ST  TIME  NODES
123456   shared     my_first_job  username  R   0:02  1
```

**Status codes:**
* `PD` = Pending (waiting in queue)
* `R` = Running  ← Your job is active!
* `CG` = Completing
* (No output) = Job finished

**Refresh every few seconds:**
```bash
watch -n 5 squeue -u $USER
```

Press `Ctrl+C` to stop watching.

### Step 7: Wait for Completion

Your job will take ~5-10 minutes. When `squeue` shows nothing, it's done!

### Step 8: Check the Results

```bash
ls -lh
```

You should see:
```text
argon.in
submit_argon.sh
job_123456.out         ← Job log
job_123456.err         ← Errors (should be empty)
thermo.out             ← Temperature/energy data
argon.lammpstrj        ← Trajectory
final_state.data       ← Final configuration
```

**Look at the job output:**
```bash
cat job_123456.out
```

You should see your echo messages plus LAMMPS output.

**Check for errors:**
```bash
cat job_123456.err
```

If this file is empty, everything worked!

---

## 4. What If It Didn't Work?

### Problem 1: "Submitted batch job" but nothing happens

**Check:** Is your job stuck in the queue?
```bash
squeue -u $USER
```

If status is `PD` (pending) for >5 minutes, the cluster might be busy. Be patient.

### Problem 2: Job finishes immediately (< 1 second)

**Check the error file:**
```bash
cat job_*.err
```

**Common errors:**

**"lmp: command not found"**
```text
FIX: Check your module load commands. Make sure you have:
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310
```

**"Cannot open input script argon.in"**
```text
FIX: LAMMPS can't find your input file.
Check: ls -lh argon.in
Make sure argon.in is in the same directory as submit_argon.sh
```

**"Invalid command" or "Unknown command"**
```text
FIX: Typo in your LAMMPS input file.
Open argon.in and check for spelling mistakes.
```

### Problem 3: Job runs but no output files

**Check:** Did the job actually finish?
```bash
tail -20 job_*.out
```

Look for "Job completed" at the end. If you don't see it, the job was killed (probably hit time limit).

---

## 5. Useful SLURM Commands

| Command | What It Does | Example |
|---------|--------------|---------|
| `sbatch script.sh` | Submit a job | `sbatch submit_argon.sh` |
| `squeue -u $USER` | Check your jobs | See what's running |
| `scancel 123456` | Cancel job 123456 | `scancel 123456` |
| `scancel -u $USER` | Cancel ALL your jobs | Use with caution! |
| `scontrol show job 123456` | Detailed job info | Why is it pending? |

---

## 6. Running Longer Jobs

For production runs (hours to days), modify the script:

```bash
nano submit_production.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=argon_production
#SBATCH --account=chm250117
#SBATCH --partition=standard          ← Use 'standard' for longer jobs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8           ← More cores for speed
#SBATCH --time=24:00:00               ← 24 hours
#SBATCH --output=prod_%j.out

module purge
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310

# Run from SCRATCH for fast I/O
cd $SCRATCH/lammps_tutorial

# Run with MPI for parallel execution
mpirun -np $SLURM_NTASKS lmp -in argon.in

echo "Production run complete"
```

---

## 7. Exercise: Run at Different Temperature

**Challenge:** Submit TWO jobs running Argon at 300 K and 700 K simultaneously.

<details>
<summary>Hint</summary>

You need:
1. Two input files: `argon_300K.in` and `argon_700K.in`
2. Two submission scripts: `submit_300K.sh` and `submit_700K.sh`

Make sure the output filenames are different in each!
</details>

<details>
<summary>Solution</summary>

```bash
# Create 300K version
cp argon.in argon_300K.in
nano argon_300K.in
# Change: velocity all create 2.50 ... and fix nvt temp 2.50 2.50 ...
# Change: log thermo_300K.out

# Create 700K version  
cp argon.in argon_700K.in
nano argon_700K.in
# Change: velocity all create 5.84 ... and fix nvt temp 5.84 5.84 ...
# Change: log thermo_700K.out

# Create submission scripts
cp submit_argon.sh submit_300K.sh
nano submit_300K.sh
# Change: lmp -in argon_300K.in

cp submit_argon.sh submit_700K.sh
nano submit_700K.sh
# Change: lmp -in argon_700K.in

# Submit both!
sbatch submit_300K.sh
sbatch submit_700K.sh

# Check both are running
squeue -u $USER
```
</details>

---

## 8. Class Examples

Pre-made scripts are available:
```bash
ls /anvil/projects/x-chm250117/class_examples/
```

Copy a complete example:
```bash
cp /anvil/projects/x-chm250117/class_examples/submit_argon.sh ./
cp /anvil/projects/x-chm250117/class_examples/argon_nvt.in ./
```

---

## 9. Summary Checklist

After this lab, you should have:

- [ ] Created a SLURM submission script
- [ ] Successfully submitted a job with `sbatch`
- [ ] Monitored job status with `squeue`
- [ ] Retrieved output files (`thermo.out`, etc.)
- [ ] Understand what each `#SBATCH` directive does

**Next Tutorial:** Analyzing LAMMPS output with Python (Tutorial 9)

---

## 10. Quick Reference Card

**Essential commands:**
```bash
# Submit job
sbatch script.sh

# Check status
squeue -u $USER

# Cancel job
scancel JOB_ID

# View output (while running)
tail -f job_*.out

# Check errors
cat job_*.err
```

**Job script template:**
```bash
#!/bin/bash
#SBATCH --account=chm250117
#SBATCH --partition=shared
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4

module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310

lmp -in input.in
```

---

## 11. Further Reading

* **[Anvil User Guide](https://www.rcac.purdue.edu/knowledge/anvil/run)** — Running jobs
* **[Anvil LAMMPS Documentation](https://www.rcac.purdue.edu/knowledge/anvil/software/installing_applications/lammps/provided_module)** — LAMMPS on Anvil
* **[SLURM Quick Start](https://slurm.schedmd.com/quickstart.html)** — Official SLURM guide
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