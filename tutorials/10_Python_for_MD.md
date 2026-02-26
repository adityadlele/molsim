# Tutorial 10: Analyzing LAMMPS Output with Python

**Objective:** Learn to process LAMMPS thermodynamic output files, detect equilibration, and calculate ensemble averages with proper uncertainty quantification.

:::{note} Prerequisites
* Tutorial 7: LAMMPS Input Files
* Tutorial 8: Running LAMMPS on Anvil
* Tutorial 4: Python and Jupyter on Anvil
* Basic understanding of statistical mechanics (mean, standard deviation)
:::

---

## 1. The Analysis Problem

You've run a LAMMPS simulation and obtained a `thermo.out` file. This file contains timestep-by-timestep snapshots of macroscopic properties (Temperature, Pressure, Energy). Now you need to answer:

1. **Has the system equilibrated?** (Is it in steady-state?)
2. **What is the average temperature?**
3. **What is the uncertainty on that average?**

These questions sound simple, but **improper analysis is the #1 source of incorrect results in MD simulations.**

### 1.1 The Two Deadly Sins

**Sin #1: Including Non-Equilibrated Data**
* Calculating averages before the system has equilibrated gives biased results
* Like measuring the temperature of water while it's still heating up

**Sin #2: Ignoring Autocorrelation**
* Consecutive data points are not independent
* Using the standard error formula $\sigma / \sqrt{N}$ without accounting for correlation severely underestimates uncertainty

This tutorial will teach you to avoid both.

---

## 2. Preparing the Simulations

Before we can analyze, we need data. We'll run **two simulations** to demonstrate the impact of equilibration.

### 2.1 Setup Directory Structure

Log in to Anvil and create a working directory:

```bash
cd $SCRATCH
mkdir lammps_analysis_tutorial
cd lammps_analysis_tutorial
```

### 2.2 Simulation 1: With Equilibration

Create `argon_equilibrated.in`:

```bash
nano argon_equilibrated.in
```

Paste this input (modified from Tutorial 7):

```lammps
# ========================================
# Argon Gas with Equilibration Phase
# Tutorial 9 - Analysis Practice
# ========================================

# ========================================
# INITIALIZATION
# ========================================
units           lj
atom_style      atomic
boundary        p p p

# ========================================
# SYSTEM DEFINITION
# ========================================
lattice         fcc 0.8
region          simbox block 0 13.57 0 13.57 0 13.57
create_box      1 simbox
create_atoms    1 box
mass            1 1.0

# Start from very low temperature (cold start)
velocity        all create 1.0 87287 dist gaussian

# ========================================
# FORCE FIELD
# ========================================
pair_style      lj/cut 3.0
pair_coeff      1 1 1.0 1.0 3.0
neighbor        0.3 bin
neigh_modify    every 1 delay 0 check yes

# ========================================
# SIMULATION SETTINGS
# ========================================
timestep        0.005
fix             1 all nvt temp 4.17 4.17 0.5

# ========================================
# OUTPUT
# ========================================
thermo_style    custom step temp press pe ke etotal density
thermo          100

# ========================================
# PHASE 1: EQUILIBRATION (20,000 steps)
# ========================================
# Heat from T*=1.0 to T*=4.17
run             20000

# ========================================
# PHASE 2: PRODUCTION (80,000 steps)
# ========================================
# Continue at equilibrated state
# Reset thermo output to new file for analysis
log             thermo_equilibrated.out
run             80000

# Save final state
write_data      final_equilibrated.data
```

### 2.3 Simulation 2: Without Equilibration

Create `argon_no_equilibration.in`:

```bash
nano argon_no_equilibration.in
```

```lammps
# ========================================
# Argon Gas WITHOUT Equilibration Phase
# Tutorial 9 - Demonstrates poor practice
# ========================================

units           lj
atom_style      atomic
boundary        p p p

lattice         fcc 0.8
region          simbox block 0 13.57 0 13.57 0 13.57
create_box      1 simbox
create_atoms    1 box
mass            1 1.0

# Start from very low temperature (cold start) - SAME AS BEFORE
velocity        all create 1.0 87287 dist gaussian

pair_style      lj/cut 3.0
pair_coeff      1 1 1.0 1.0 3.0
neighbor        0.3 bin
neigh_modify    every 1 delay 0 check yes

timestep        0.005
fix             1 all nvt temp 4.17 4.17 0.5

thermo_style    custom step temp press pe ke etotal density
thermo          100
log             thermo_no_equilibration.out

# Run production immediately (no equilibration!)
run             80000

write_data      final_no_equilibration.data
```

### 2.4 Submit Both Simulations

Create a submission script:

```bash
nano submit_analysis_jobs.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=analysis_prep
#SBATCH --account=chm250117
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=01:00:00
#SBATCH --output=prep_%j.out

module purge
module load gcc/11.2.0 openmpi/4.0.6
module load lammps/20210310

echo "Running equilibrated simulation..."
lmp -in argon_equilibrated.in

echo "Running non-equilibrated simulation..."
lmp -in argon_no_equilibration.in

echo "Both simulations complete. Ready for analysis."
```

Submit:
```bash
sbatch submit_analysis_jobs.sh
```

Wait for completion (~10-20 minutes). Check with:
```bash
squeue -u $USER
```

Once done, verify you have the output files:
```bash
ls -lh thermo*.out
```

---

## 3. Launching Jupyter for Analysis

Now we switch to Python for analysis.

### 3.1 Start Jupyter Session

1. Go to [Anvil OnDemand](https://ondemand.anvil.rcac.purdue.edu)
2. Click **Interactive Apps** → **Jupyter Notebook**
3. Settings:
   * Account: `chm250117`
   * Partition: `shared`
   * Number of hours: `2`
4. Launch and connect

### 3.2 Navigate to Working Directory

In Jupyter, click **New** → **Terminal**, then:
```bash
cd $SCRATCH/lammps_analysis_tutorial
```

Close the terminal. In Jupyter file browser, navigate to `lammps_analysis_tutorial`.

### 3.3 Create Analysis Notebook

Click **New** → **Python [conda env: molsimclass]**

Rename the notebook to `LAMMPS_Analysis.ipynb`.

:::{tip} Using Class Examples
Sample LAMMPS output files are available in the class directory:
```bash
ls /anvil/projects/x-chm250117/class_examples/
```
You can copy example `thermo.out` files to practice analysis without running your own simulations first.
:::

---

## 4. Loading and Parsing LAMMPS Output

### 4.1 Understanding `thermo.out` Format

Open `thermo_equilibrated.out` in a text editor to see the structure:

```text
LAMMPS (version)
...
Step Temp Press PotEng KinEng TotEng Density 
0 1.0000000 4.0084736 -2.6604462 1.4977500 -1.1626962 0.8000000 
100 1.5234567 6.1234567 -2.5432198 2.2834567 -0.2597631 0.8000000 
200 2.0123456 8.2345678 -2.4321987 3.0123456 0.5801469 0.8000000 
...
```

The file has:
* Header lines (starting with `LAMMPS`, `#`, etc.)
* Column names (`Step Temp Press ...`)
* Data rows (numerical values)

### 4.2 Python: Loading the Data

**Cell 1: Import Libraries**

```python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

# Configure plotting
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.figsize'] = (12, 5)
plt.rcParams['font.size'] = 11

print("Libraries loaded successfully")
```

**Cell 2: Parse LAMMPS Output**

```python
def load_lammps_thermo(filename):
    """
    Parse LAMMPS thermo output file.
    
    Parameters
    ----------
    filename : str
        Path to thermo.out file
    
    Returns
    -------
    pandas.DataFrame
        Thermodynamic data with columns as headers
    """
    # Read file, skip header lines, use first data row as column names
    data = pd.read_csv(
        filename,
        sep=r'\s+',           # Whitespace-separated
        comment='#',          # Ignore lines starting with #
        skiprows=lambda x: x < 2  # Skip LAMMPS header
    )
    
    # First row might be column names, check and reset if needed
    if data.iloc[0, 0] == 'Step':
        data.columns = data.iloc[0]
        data = data[1:].reset_index(drop=True)
        data = data.apply(pd.to_numeric, errors='coerce')
    
    return data

# Test loading
data_eq = load_lammps_thermo('thermo_equilibrated.out')
data_no_eq = load_lammps_thermo('thermo_no_equilibration.out')

print("Equilibrated data shape:", data_eq.shape)
print("Non-equilibrated data shape:", data_no_eq.shape)
print("\nColumn names:", list(data_eq.columns))
print("\nFirst few rows:")
print(data_eq.head())
```

---

## 5. Visual Equilibration Detection

The first step in analysis is **always** to plot the raw data.

### 5.1 Temperature Time Series

**Cell 3: Plot Temperature Evolution**

```python
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: With Equilibration
axes[0].plot(data_eq['Step'], data_eq['Temp'], 
             color='steelblue', alpha=0.7, lw=0.8)
axes[0].axhline(4.17, color='red', linestyle='--', 
                lw=2, label='Target T* = 4.17')
axes[0].set_xlabel('Step')
axes[0].set_ylabel('Temperature (T*)')
axes[0].set_title('With Equilibration Phase')
axes[0].legend()
axes[0].grid(alpha=0.3)

# Plot 2: Without Equilibration
axes[1].plot(data_no_eq['Step'], data_no_eq['Temp'], 
             color='darkorange', alpha=0.7, lw=0.8)
axes[1].axhline(4.17, color='red', linestyle='--', 
                lw=2, label='Target T* = 4.17')
axes[1].set_xlabel('Step')
axes[1].set_ylabel('Temperature (T*)')
axes[1].set_title('WITHOUT Equilibration Phase')
axes[1].legend()
axes[1].grid(alpha=0.3)

plt.suptitle('Temperature Evolution Comparison', 
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.show()

print("\nObservation:")
print("  Left plot: System starts equilibrated (from previous run)")
print("  Right plot: System starts cold (T*=1.0) and heats up")
print("  → Must discard initial transient data!")
```

### 5.2 Running Average (Moving Average)

A **running average** smooths out fluctuations to reveal the underlying trend.

**Cell 4: Compute and Plot Running Averages**

```python
def running_average(data, window):
    """
    Compute running (moving) average using convolution.
    
    Parameters
    ----------
    data : array-like
        Input time series
    window : int
        Window size (number of points to average)
    
    Returns
    -------
    ndarray
        Smoothed data (shorter by window-1 points)
    """
    kernel = np.ones(window) / window
    return np.convolve(data, kernel, mode='valid')

# Calculate running averages
window = 200
temp_eq_smooth = running_average(data_eq['Temp'].values, window)
temp_no_eq_smooth = running_average(data_no_eq['Temp'].values, window)

# Adjust x-axis (convolution shortens array)
steps_eq_smooth = data_eq['Step'].values[window-1:]
steps_no_eq_smooth = data_no_eq['Step'].values[window-1:]

# Plot
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: With Equilibration
axes[0].plot(data_eq['Step'], data_eq['Temp'], 
             color='steelblue', alpha=0.3, lw=0.5, label='Raw')
axes[0].plot(steps_eq_smooth, temp_eq_smooth, 
             color='navy', lw=2.5, label=f'Running Avg (window={window})')
axes[0].axhline(4.17, color='red', linestyle='--', lw=1.5, alpha=0.7)
axes[0].set_xlabel('Step')
axes[0].set_ylabel('Temperature (T*)')
axes[0].set_title('With Equilibration')
axes[0].legend()
axes[0].set_ylim(0, 6)

# Plot 2: Without Equilibration
axes[1].plot(data_no_eq['Step'], data_no_eq['Temp'], 
             color='darkorange', alpha=0.3, lw=0.5, label='Raw')
axes[1].plot(steps_no_eq_smooth, temp_no_eq_smooth, 
             color='darkred', lw=2.5, label=f'Running Avg (window={window})')
axes[1].axhline(4.17, color='red', linestyle='--', lw=1.5, alpha=0.7)
axes[1].set_xlabel('Step')
axes[1].set_ylabel('Temperature (T*)')
axes[1].set_title('WITHOUT Equilibration')
axes[1].legend()
axes[1].set_ylim(0, 6)

plt.suptitle('Running Average Reveals Equilibration Drift', 
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.show()

print("\nInterpretation:")
print("  If running average is horizontal → System is equilibrated")
print("  If running average slopes up/down → Still equilibrating (transient phase)")
```

---

## 6. Quantitative Equilibration Detection

Visual inspection works for obvious cases, but we need a **quantitative criterion** for automation.

### 6.1 Block Analysis Method

**Concept:** Divide the simulation into blocks. If the mean of each block is roughly constant, the system is equilibrated.

**Cell 5: Block Analysis for Equilibration**

```python
def detect_equilibration(data, n_blocks=5):
    """
    Detect equilibration using block analysis.
    
    Strategy: Divide data into blocks. If std(block_means) is small
    relative to overall std, system is equilibrated.
    
    Parameters
    ----------
    data : array-like
        Time series (e.g., temperature)
    n_blocks : int
        Number of blocks to divide data into
    
    Returns
    -------
    dict
        Contains block means, overall stats, and equilibration assessment
    """
    n_total = len(data)
    block_size = n_total // n_blocks
    
    block_means = []
    for i in range(n_blocks):
        start = i * block_size
        end = (i + 1) * block_size if i < n_blocks - 1 else n_total
        block_means.append(np.mean(data[start:end]))
    
    overall_mean = np.mean(data)
    overall_std = np.std(data, ddof=1)
    block_std = np.std(block_means, ddof=1)
    
    # Criterion: if block_std < 0.1 * overall_std, likely equilibrated
    is_equilibrated = (block_std < 0.1 * overall_std)
    
    return {
        'block_means': block_means,
        'overall_mean': overall_mean,
        'overall_std': overall_std,
        'block_std': block_std,
        'equilibrated': is_equilibrated
    }

# Test on both datasets
result_eq = detect_equilibration(data_eq['Temp'].values)
result_no_eq = detect_equilibration(data_no_eq['Temp'].values)

print("="*60)
print("EQUILIBRATION DETECTION RESULTS")
print("="*60)
print("\nWith Equilibration Phase:")
print(f"  Block std: {result_eq['block_std']:.4f}")
print(f"  Overall std: {result_eq['overall_std']:.4f}")
print(f"  Ratio: {result_eq['block_std']/result_eq['overall_std']:.3f}")
print(f"  Equilibrated: {result_eq['equilibrated']}")

print("\nWithout Equilibration Phase:")
print(f"  Block std: {result_no_eq['block_std']:.4f}")
print(f"  Overall std: {result_no_eq['overall_std']:.4f}")
print(f"  Ratio: {result_no_eq['block_std']/result_no_eq['overall_std']:.3f}")
print(f"  Equilibrated: {result_no_eq['equilibrated']}")
print("="*60)
```

### 6.2 Manual Truncation

For the non-equilibrated data, we need to **discard the initial transient**. From the plot, we can see the system reaches steady state around step 20,000.

**Cell 6: Truncate Non-Equilibrated Data**

```python
# Visually identified equilibration point: step 20000
# Find corresponding index
equilibration_cutoff = 20000
cutoff_index = data_no_eq[data_no_eq['Step'] >= equilibration_cutoff].index[0]

print(f"Discarding first {cutoff_index} rows (steps 0-{equilibration_cutoff})")

# Create truncated dataset
data_no_eq_truncated = data_no_eq.iloc[cutoff_index:].reset_index(drop=True)

print(f"Original length: {len(data_no_eq)} rows")
print(f"Truncated length: {len(data_no_eq_truncated)} rows")
print(f"Discarded: {100 * cutoff_index / len(data_no_eq):.1f}% of data")

# Verify equilibration on truncated data
result_truncated = detect_equilibration(data_no_eq_truncated['Temp'].values)
print(f"\nAfter truncation:")
print(f"  Equilibrated: {result_truncated['equilibrated']}")
```

---

## 7. Calculating Ensemble Averages

Now that we have equilibrated data, we can calculate meaningful averages.

**Cell 7: Compute Average Properties**

```python
def calculate_ensemble_averages(data):
    """
    Calculate mean and standard deviation for all properties.
    
    Parameters
    ----------
    data : pandas.DataFrame
        Thermodynamic data
    
    Returns
    -------
    pandas.DataFrame
        Summary statistics
    """
    properties = ['Temp', 'Press', 'PotEng', 'KinEng', 'TotEng', 'Density']
    results = []
    
    for prop in properties:
        if prop in data.columns:
            mean_val = np.mean(data[prop])
            std_val = np.std(data[prop], ddof=1)
            results.append({
                'Property': prop,
                'Mean': mean_val,
                'Std': std_val,
                'Std/Mean (%)': 100 * std_val / abs(mean_val) if mean_val != 0 else np.nan
            })
    
    return pd.DataFrame(results)

# Calculate for both datasets
summary_eq = calculate_ensemble_averages(data_eq)
summary_no_eq_truncated = calculate_ensemble_averages(data_no_eq_truncated)

print("="*70)
print("ENSEMBLE AVERAGES: With Equilibration Phase")
print("="*70)
print(summary_eq.to_string(index=False))

print("\n" + "="*70)
print("ENSEMBLE AVERAGES: Without Equilibration (After Truncation)")
print("="*70)
print(summary_no_eq_truncated.to_string(index=False))
```

---

## 8. Uncertainty Quantification: The Autocorrelation Problem

### 8.1 Why Standard Error is Wrong

The naive standard error formula is:
$$\text{SE} = \frac{\sigma}{\sqrt{N}}$$

This assumes **independent** samples. But consecutive MD snapshots are highly correlated!

### 8.2 Statistical Inefficiency

We define the **statistical inefficiency** $s$ as:
$$s = 1 + 2 \sum_{t=1}^{\infty} C(t)$$

where $C(t)$ is the autocorrelation function. The **effective sample size** is:
$$N_{\text{eff}} = \frac{N}{s}$$

The **correct standard error** is:
$$\text{SE}_{\text{corrected}} = \frac{\sigma}{\sqrt{N_{\text{eff}}}} = \frac{\sigma}{\sqrt{N/s}} = \frac{\sigma \sqrt{s}}{\sqrt{N}}$$

### 8.3 Block Averaging Method

**The robust practical method:** Divide data into $k$ blocks, compute block means, then:
$$\text{SE} = \frac{\sigma_{\text{blocks}}}{\sqrt{k}}$$

where $\sigma_{\text{blocks}}$ is the standard deviation of the block means.

**Cell 8: Block Averaging for Uncertainty**

```python
def block_average_uncertainty(data, n_blocks=10):
    """
    Calculate standard error using block averaging.
    
    This accounts for autocorrelation by treating block means
    as independent samples.
    
    Parameters
    ----------
    data : array-like
        Time series
    n_blocks : int
        Number of blocks
    
    Returns
    -------
    dict
        Mean, standard error, and confidence interval
    """
    n_total = len(data)
    block_size = n_total // n_blocks
    
    block_means = []
    for i in range(n_blocks):
        start = i * block_size
        end = (i + 1) * block_size if i < n_blocks - 1 else n_total
        block_means.append(np.mean(data[start:end]))
    
    overall_mean = np.mean(block_means)
    block_std = np.std(block_means, ddof=1)
    standard_error = block_std / np.sqrt(n_blocks)
    
    # 95% confidence interval (t-distribution for small samples)
    t_critical = stats.t.ppf(0.975, df=n_blocks-1)
    ci_95 = t_critical * standard_error
    
    # Naive SE for comparison
    naive_se = np.std(data, ddof=1) / np.sqrt(n_total)
    
    return {
        'mean': overall_mean,
        'se': standard_error,
        'ci_95': ci_95,
        'naive_se': naive_se,
        'inflation_factor': standard_error / naive_se
    }

# Calculate for temperature
temp_uncertainty_eq = block_average_uncertainty(data_eq['Temp'].values, n_blocks=10)
temp_uncertainty_no_eq = block_average_uncertainty(
    data_no_eq_truncated['Temp'].values, n_blocks=10
)

print("="*70)
print("UNCERTAINTY QUANTIFICATION: Temperature")
print("="*70)

print("\nWith Equilibration:")
print(f"  Mean:              {temp_uncertainty_eq['mean']:.4f}")
print(f"  Standard Error:    {temp_uncertainty_eq['se']:.4f}")
print(f"  95% CI:            ±{temp_uncertainty_eq['ci_95']:.4f}")
print(f"  Naive SE:          {temp_uncertainty_eq['naive_se']:.4f}")
print(f"  Inflation Factor:  {temp_uncertainty_eq['inflation_factor']:.2f}x")
print(f"  → T = {temp_uncertainty_eq['mean']:.4f} ± {temp_uncertainty_eq['ci_95']:.4f}")

print("\nWithout Equilibration (Truncated):")
print(f"  Mean:              {temp_uncertainty_no_eq['mean']:.4f}")
print(f"  Standard Error:    {temp_uncertainty_no_eq['se']:.4f}")
print(f"  95% CI:            ±{temp_uncertainty_no_eq['ci_95']:.4f}")
print(f"  Naive SE:          {temp_uncertainty_no_eq['naive_se']:.4f}")
print(f"  Inflation Factor:  {temp_uncertainty_no_eq['inflation_factor']:.2f}x")
print(f"  → T = {temp_uncertainty_no_eq['mean']:.4f} ± {temp_uncertainty_no_eq['ci_95']:.4f}")

print("\n" + "="*70)
print("INTERPRETATION:")
print("  The 'Inflation Factor' shows how much the naive formula")
print("  underestimates uncertainty due to autocorrelation.")
print("  Values of 2-5x are typical for MD simulations.")
print("="*70)
```

---

## 9. Complete Analysis Function

**Cell 9: Automated Analysis Pipeline**

```python
def analyze_lammps_run(filename, equilibration_fraction=0.2, n_blocks=10):
    """
    Complete analysis pipeline for LAMMPS thermo output.
    
    Parameters
    ----------
    filename : str
        Path to thermo.out file
    equilibration_fraction : float
        Fraction of data to discard as equilibration (default: 0.2 = 20%)
    n_blocks : int
        Number of blocks for uncertainty estimation
    
    Returns
    -------
    dict
        Contains data, summary statistics, and uncertainty estimates
    """
    # Load data
    data = load_lammps_thermo(filename)
    
    # Discard equilibration
    cutoff = int(len(data) * equilibration_fraction)
    data_production = data.iloc[cutoff:].reset_index(drop=True)
    
    print(f"Loaded: {filename}")
    print(f"Total steps: {len(data)}")
    print(f"Discarded: {cutoff} steps ({equilibration_fraction*100:.0f}%)")
    print(f"Production: {len(data_production)} steps")
    
    # Calculate properties
    properties = ['Temp', 'Press', 'PotEng', 'KinEng', 'TotEng']
    results = {}
    
    for prop in properties:
        if prop in data_production.columns:
            uncertainty = block_average_uncertainty(
                data_production[prop].values, n_blocks
            )
            results[prop] = uncertainty
    
    return {
        'data_full': data,
        'data_production': data_production,
        'properties': results,
        'n_production': len(data_production)
    }

# Run complete analysis
analysis_eq = analyze_lammps_run('thermo_equilibrated.out', 
                                  equilibration_fraction=0.0, 
                                  n_blocks=10)
analysis_no_eq = analyze_lammps_run('thermo_no_equilibration.out', 
                                      equilibration_fraction=0.25, 
                                      n_blocks=10)

print("\n" + "="*70)
print("FINAL RESULTS SUMMARY")
print("="*70)

print("\nTemperature (T*):")
print(f"  Equilibrated:     {analysis_eq['properties']['Temp']['mean']:.4f} "
      f"± {analysis_eq['properties']['Temp']['ci_95']:.4f}")
print(f"  Non-equilibrated: {analysis_no_eq['properties']['Temp']['mean']:.4f} "
      f"± {analysis_no_eq['properties']['Temp']['ci_95']:.4f}")

print("\nPressure (P*):")
print(f"  Equilibrated:     {analysis_eq['properties']['Press']['mean']:.4f} "
      f"± {analysis_eq['properties']['Press']['ci_95']:.4f}")
print(f"  Non-equilibrated: {analysis_no_eq['properties']['Press']['mean']:.4f} "
      f"± {analysis_no_eq['properties']['Press']['ci_95']:.4f}")
```

---

## 10. Visualizing the Final Results

**Cell 10: Create Publication-Quality Summary Plot**

```python
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Define properties and their labels
properties = [
    ('Temp', 'Temperature (T*)', 4.17),
    ('Press', 'Pressure (P*)', None),
    ('TotEng', 'Total Energy (E*)', None),
    ('Density', 'Density (ρ*)', 0.8)
]

for idx, (prop, label, target) in enumerate(properties):
    row, col = idx // 2, idx % 2
    ax = axes[row, col]
    
    # Plot raw data
    if prop in analysis_eq['data_production'].columns:
        ax.plot(analysis_eq['data_production']['Step'], 
                analysis_eq['data_production'][prop],
                color='steelblue', alpha=0.5, lw=0.6, label='Equilibrated')
    
    if prop in analysis_no_eq['data_production'].columns:
        ax.plot(analysis_no_eq['data_production']['Step'], 
                analysis_no_eq['data_production'][prop],
                color='darkorange', alpha=0.5, lw=0.6, label='Non-equilibrated')
    
    # Add mean lines
    if prop in analysis_eq['properties']:
        mean_eq = analysis_eq['properties'][prop]['mean']
        ci_eq = analysis_eq['properties'][prop]['ci_95']
        ax.axhline(mean_eq, color='navy', lw=2, ls='-', 
                   label=f'Mean ± CI: {mean_eq:.3f} ± {ci_eq:.3f}')
        ax.axhspan(mean_eq - ci_eq, mean_eq + ci_eq, 
                   alpha=0.2, color='steelblue')
    
    # Add target line if applicable
    if target is not None:
        ax.axhline(target, color='red', lw=1.5, ls='--', 
                   alpha=0.7, label=f'Target: {target}')
    
    ax.set_xlabel('Step')
    ax.set_ylabel(label)
    ax.set_title(label)
    ax.legend(fontsize=9, loc='best')
    ax.grid(alpha=0.3)

plt.suptitle('LAMMPS Output Analysis: Equilibration Comparison', 
             fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig('lammps_analysis_summary.png', dpi=150, bbox_inches='tight')
plt.show()

print("Summary plot saved as: lammps_analysis_summary.png")
```

---

## 11. Exercise: Effect of Simulation Length

**Task:** Investigate how simulation length affects the uncertainty on estimated properties.

**Cell 11: Simulation Length Study**

```python
def uncertainty_vs_length(data, property_name, max_blocks=20):
    """
    Calculate how uncertainty decreases with simulation length.
    
    Parameters
    ----------
    data : array-like
        Property time series
    property_name : str
        Name for plotting
    max_blocks : int
        Maximum number of blocks to test
    
    Returns
    -------
    dict
        Contains arrays of n_blocks, standard errors, and confidence intervals
    """
    n_blocks_array = np.arange(2, max_blocks + 1)
    se_array = []
    ci_array = []
    
    for n_blocks in n_blocks_array:
        result = block_average_uncertainty(data, n_blocks)
        se_array.append(result['se'])
        ci_array.append(result['ci_95'])
    
    return {
        'n_blocks': n_blocks_array,
        'se': np.array(se_array),
        'ci': np.array(ci_array)
    }

# Calculate for temperature
temp_data = analysis_eq['data_production']['Temp'].values
length_study = uncertainty_vs_length(temp_data, 'Temperature', max_blocks=20)

# Plot
fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(length_study['n_blocks'], length_study['se'], 
        'o-', color='steelblue', lw=2, markersize=6, label='Standard Error')
ax.plot(length_study['n_blocks'], length_study['ci'], 
        's-', color='darkorange', lw=2, markersize=6, label='95% CI')

# Add theoretical 1/sqrt(N) line
n_ref = length_study['n_blocks'][0]
se_ref = length_study['se'][0]
theoretical = se_ref * np.sqrt(n_ref / length_study['n_blocks'])
ax.plot(length_study['n_blocks'], theoretical, 
        '--', color='gray', lw=2, alpha=0.7, label='Theoretical (1/√N)')

ax.set_xlabel('Number of Blocks (proxy for simulation length)', fontsize=12)
ax.set_ylabel('Uncertainty', fontsize=12)
ax.set_title('How Uncertainty Decreases with Simulation Length', fontsize=14)
ax.legend(fontsize=11)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.show()

print("Interpretation:")
print("  Uncertainty decreases as ~1/√N with simulation length")
print("  To halve uncertainty, you need 4x longer simulation")
print("  To reduce by 10x, you need 100x longer simulation")
```

---

## 12. Best Practices Summary

**DO:**
- ✅ Always plot raw data before calculating averages
- ✅ Use running averages to visually detect equilibration
- ✅ Discard equilibration phase (typically 10-25% of run)
- ✅ Use block averaging for uncertainty quantification
- ✅ Report results as: Mean ± 95% CI
- ✅ Run multiple independent simulations for robustness

**DON'T:**
- ❌ Calculate averages before equilibration
- ❌ Use naive standard error $\sigma/\sqrt{N}$ without checking autocorrelation
- ❌ Trust a single simulation run
- ❌ Ignore large fluctuations or drift in running average
- ❌ Mix equilibration and production data

---

## 13. Exporting Results

**Cell 12: Save Results to CSV**

```python
# Create summary DataFrame
summary_data = []

for prop in ['Temp', 'Press', 'PotEng', 'TotEng']:
    if prop in analysis_eq['properties']:
        eq_result = analysis_eq['properties'][prop]
        no_eq_result = analysis_no_eq['properties'][prop]
        
        summary_data.append({
            'Property': prop,
            'With_Eq_Mean': eq_result['mean'],
            'With_Eq_CI': eq_result['ci_95'],
            'No_Eq_Mean': no_eq_result['mean'],
            'No_Eq_CI': no_eq_result['ci_95'],
            'Difference': abs(eq_result['mean'] - no_eq_result['mean'])
        })

summary_df = pd.DataFrame(summary_data)

# Save to CSV
summary_df.to_csv('lammps_analysis_results.csv', index=False)

print("Results saved to: lammps_analysis_results.csv")
print("\nSummary:")
print(summary_df.to_string(index=False))
```

---

## 14. Further Reading

* **Allen & Tildesley, "Computer Simulation of Liquids"** — Chapter 6: Statistical Mechanics
* **Frenkel & Smit, "Understanding Molecular Simulation"** — Chapter 4: Molecular Dynamics Simulations
* **Flyvbjerg & Petersen (1989)** — "Error estimates on averages of correlated data" (classic paper on block averaging)
* **[LAMMPS Howto: Output](https://docs.lammps.org/Howto_output.html)** — Official guide to LAMMPS output formats

:::{tip} Automating Analysis
The code in this tutorial can be packaged into a Python script and run from the command line:
```bash
python analyze_lammps.py thermo.out --equilibration 0.2 --blocks 10
```
This is useful for processing many simulations in batch.
:::

:::{warning} When Equilibration Fails
If your running average never flattens out (continues drifting after 100,000+ steps), your system may have:
* Multiple metastable states (e.g., crystal vs. liquid)
* Very slow relaxation (glass-forming systems)
* Incorrect force field parameters
* Box size too small (finite-size effects)

In these cases, advanced techniques (enhanced sampling, longer simulations, or different ensembles) are needed.
:::

---

## 15. Next Steps

Now that you can analyze LAMMPS output:
* **Midterm Exam Preparation:** Practice running simulations at different conditions and analyzing trends
* **Advanced Analysis:** Learn to calculate radial distribution functions (RDF), mean squared displacement (MSD), and transport properties
* **Visualization:** Use tools like OVITO or VMD to visualize trajectories

This tutorial completes the core LAMMPS workflow: **Build** (Tutorial 7) → **Run** (Tutorial 8) → **Analyze** (Tutorial 9).