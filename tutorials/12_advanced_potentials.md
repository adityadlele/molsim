# Tutorial 12: Advanced MD Potentials in LAMMPS

**Objective:** Learn to use professional force field parameter files, simulate metallic systems with EAM, and run reactive chemistry simulations with ReaxFF.

---

## Using Pre-Run Simulations

For this tutorial, we have pre-run all simulations in a shared class directory. You can either:

1. **Use pre-run results** (recommended for in-class work): Copy outputs and run analysis
2. **Run simulations yourself** (for practice later): Use the LAMMPS input files and geometry builders provided

**Pre-run simulation location:**
```bash
/anvil/projects/x-chm250117/class_examples/tutorial12_completed/
```

**Directory structure:**
```
/anvil/projects/x-chm250117/class_examples/tutorial12_completed/
├── part1_opls_methanol/
│   ├── methanol.data
│   ├── methanol.in
│   ├── methanol.log
│   ├── methanol.lammpstrj
│   ├── methanol_rdf.dat
│   └── methanol_final.data
├── part2_eam_copper/
│   ├── copper_melt.in
│   ├── copper_melt.log
│   ├── copper_heat.lammpstrj
│   ├── msd.dat
│   └── copper_final.data
└── part3_reaxff_combustion/
    ├── combustion.data
    ├── combustion.in
    ├── combustion.log
    ├── combustion.lammpstrj
    └── species.out
```

---

## File Organization

For your own work, create a well-organized directory structure on Anvil:

```bash
cd $SCRATCH
mkdir tutorial12
cd tutorial12

# Create subdirectories for each part
mkdir part1_opls_methanol
mkdir part2_eam_copper
mkdir part3_reaxff_combustion
```

---

## Part 1: Using Force Field Parameter Files (OPLS)

### 1.1 The Challenge: Real Molecules

In Tutorial 11, we manually specified all force field parameters in the LAMMPS input file. For real molecules (drugs, proteins, polymers), this is impractical:
- Hundreds of atom types
- Thousands of bond/angle/dihedral parameters
- Risk of errors and inconsistencies

**Solution:** Use pre-validated force field parameter files.

### 1.2 System: Liquid Methanol

Methanol (CH₃OH) is a simple organic molecule that requires:
- Multiple atom types (C, H on C, O, H on O)
- Bonds, angles, dihedrals
- Partial charges
- LJ parameters

We'll use **OPLS-AA** (All-Atom) force field.

### 1.3 Building the System

The data file `methanol.data` is already provided in the shared directory:

```bash
ls /anvil/projects/x-chm250117/class_examples/tutorial12_completed/part1_opls_methanol/methanol.data
```

:::{dropdown} Python Script to Build Methanol System (for reference)
```python
import numpy as np

# Parameters
n_molecules = 125        # 5x5x5 grid
box_size = 24.0          # Angstroms (density ~0.8 g/cm³)

# Methanol geometry (CH3OH)
# Atom positions in molecule-fixed frame
# C at origin, O along +x
c_pos = np.array([0.0, 0.0, 0.0])
o_pos = np.array([1.43, 0.0, 0.0])
h_oh = o_pos + np.array([0.945, 0.0, 0.0])  # OH hydrogen

# CH3 hydrogens (tetrahedral)
theta = np.radians(109.47)
h1_ch3 = c_pos + 1.09 * np.array([-np.cos(theta), np.sin(theta), 0.0])
h2_ch3 = c_pos + 1.09 * np.array([-np.cos(theta), -np.sin(theta)*0.5, np.sqrt(3)/2*np.sin(theta)])
h3_ch3 = c_pos + 1.09 * np.array([-np.cos(theta), -np.sin(theta)*0.5, -np.sqrt(3)/2*np.sin(theta)])

print("Building methanol liquid box...")
print(f"  Molecules: {n_molecules}")
print(f"  Box size: {box_size} Å")
print(f"  Atoms per molecule: 6 (C, O, 4H)")

with open('methanol.data', 'w') as f:
    f.write('OPLS Methanol System\n\n')
    f.write(f'{n_molecules*6} atoms\n')
    f.write(f'{n_molecules*5} bonds\n')
    f.write(f'{n_molecules*7} angles\n')
    f.write(f'{n_molecules*3} dihedrals\n\n')
    
    f.write('4 atom types  # 1=C(CH3), 2=O, 3=H(CH3), 4=H(OH)\n')
    f.write('3 bond types  # 1=C-O, 2=C-H, 3=O-H\n')
    f.write('3 angle types\n')
    f.write('1 dihedral types\n\n')
    
    f.write(f'0.0 {box_size} xlo xhi\n')
    f.write(f'0.0 {box_size} ylo yhi\n')
    f.write(f'0.0 {box_size} zlo zhi\n\n')
    
    f.write('Masses\n\n')
    f.write('1 12.011  # C\n')
    f.write('2 15.999  # O\n')
    f.write('3 1.008   # H(CH3)\n')
    f.write('4 1.008   # H(OH)\n\n')
    
    # OPLS-AA partial charges for methanol
    q_c = 0.145
    q_o = -0.683
    q_h_ch3 = 0.040
    q_h_oh = 0.418
    
    f.write('Atoms  # full\n\n')
    
    # Place molecules on grid with random orientations
    nx = int(round(n_molecules**(1/3)))
    spacing = box_size / nx
    
    atom_id = 1
    mol_id = 1
    
    for i in range(nx):
        for j in range(nx):
            for k in range(nx):
                # Random rotation
                phi = np.random.rand() * 2 * np.pi
                theta_r = np.random.rand() * np.pi
                
                cos_t, sin_t = np.cos(theta_r), np.sin(theta_r)
                cos_p, sin_p = np.cos(phi), np.sin(phi)
                
                # Simple rotation matrix
                rot = np.array([
                    [cos_p*cos_t, -sin_p, cos_p*sin_t],
                    [sin_p*cos_t, cos_p, sin_p*sin_t],
                    [-sin_t, 0, cos_t]
                ])
                
                # Molecule center
                center = np.array([
                    (i + 0.5) * spacing + np.random.uniform(-0.2, 0.2),
                    (j + 0.5) * spacing + np.random.uniform(-0.2, 0.2),
                    (k + 0.5) * spacing + np.random.uniform(-0.2, 0.2)
                ])
                
                # Rotate and place atoms
                atoms_local = [c_pos, o_pos, h1_ch3, h2_ch3, h3_ch3, h_oh]
                types = [1, 2, 3, 3, 3, 4]
                charges = [q_c, q_o, q_h_ch3, q_h_ch3, q_h_ch3, q_h_oh]
                
                for pos, atype, q in zip(atoms_local, types, charges):
                    rotated_pos = center + rot @ pos
                    f.write(f'{atom_id} {mol_id} {atype} {q} {rotated_pos[0]:.6f} '
                           f'{rotated_pos[1]:.6f} {rotated_pos[2]:.6f}\n')
                    atom_id += 1
                
                mol_id += 1
    
    # Write bonds
    f.write('\nBonds\n\n')
    bond_id = 1
    for mol in range(1, n_molecules + 1):
        c_atom = (mol - 1) * 6 + 1
        o_atom = c_atom + 1
        h1 = c_atom + 2
        h2 = c_atom + 3
        h3 = c_atom + 4
        h_o = c_atom + 5
        
        bonds = [
            (1, c_atom, o_atom),   # C-O
            (2, c_atom, h1),       # C-H
            (2, c_atom, h2),
            (2, c_atom, h3),
            (3, o_atom, h_o),      # O-H
        ]
        for btype, a1, a2 in bonds:
            f.write(f'{bond_id} {btype} {a1} {a2}\n')
            bond_id += 1
    
    # Write angles
    f.write('\nAngles\n\n')
    angle_id = 1
    for mol in range(1, n_molecules + 1):
        c_atom = (mol - 1) * 6 + 1
        o_atom = c_atom + 1
        h1 = c_atom + 2
        h2 = c_atom + 3
        h3 = c_atom + 4
        h_o = c_atom + 5
        
        angles = [
            (1, h1, c_atom, o_atom),  # H-C-O
            (1, h2, c_atom, o_atom),
            (1, h3, c_atom, o_atom),
            (2, c_atom, o_atom, h_o), # C-O-H
            (3, h1, c_atom, h2),      # H-C-H
            (3, h1, c_atom, h3),
            (3, h2, c_atom, h3),
        ]
        for atype, a1, a2, a3 in angles:
            f.write(f'{angle_id} {atype} {a1} {a2} {a3}\n')
            angle_id += 1
    
    # Write dihedrals
    f.write('\nDihedrals\n\n')
    dihedral_id = 1
    for mol in range(1, n_molecules + 1):
        c_atom = (mol - 1) * 6 + 1
        o_atom = c_atom + 1
        h1 = c_atom + 2
        h2 = c_atom + 3
        h3 = c_atom + 4
        h_o = c_atom + 5
        
        dihedrals = [
            (1, h1, c_atom, o_atom, h_o),  # H-C-O-H
            (1, h2, c_atom, o_atom, h_o),
            (1, h3, c_atom, o_atom, h_o),
        ]
        for dtype, a1, a2, a3, a4 in dihedrals:
            f.write(f'{dihedral_id} {dtype} {a1} {a2} {a3} {a4}\n')
            dihedral_id += 1

print(f"\n✓ Created methanol.data")
print(f"✓ Total atoms: {n_molecules * 6}")
print(f"✓ Total molecules: {n_molecules}")

import os
if os.path.exists('methanol.data'):
    size = os.path.getsize('methanol.data') / 1024
    print(f"✓ File size: {size:.1f} KB")
else:
    print("✗ ERROR: File not created!")
```
:::

### 1.4 LAMMPS Input with OPLS Parameters

The LAMMPS input file `methanol.in`:

```lammps
# Liquid Methanol - OPLS-AA Force Field
# Tutorial 12 - Part 1

units           real
atom_style      full
boundary        p p p

read_data       methanol.data

# ========================================
# OPLS-AA Force Field Parameters
# ========================================

# Non-bonded: LJ + Coulomb
pair_style      lj/cut/coul/long 10.0
# OPLS-AA parameters for methanol
pair_coeff      1 1 0.066 3.50    # C(CH3)
pair_coeff      2 2 0.170 3.12    # O
pair_coeff      3 3 0.030 2.50    # H(CH3)
pair_coeff      4 4 0.000 0.00    # H(OH) - no LJ

kspace_style    pppm 1.0e-4

# Bonds (OPLS-AA)
bond_style      harmonic
bond_coeff      1 320.0 1.410    # C-O
bond_coeff      2 340.0 1.090    # C-H
bond_coeff      3 553.0 0.945    # O-H

# Angles (OPLS-AA)
angle_style     harmonic
angle_coeff     1 50.0 109.5     # H-C-O
angle_coeff     2 55.0 108.5     # C-O-H
angle_coeff     3 33.0 107.8     # H-C-H

# Dihedrals (OPLS style)
dihedral_style  opls
dihedral_coeff  1 0.0 0.0 0.468 0.0    # H-C-O-H

# Exclusions (OPLS standard)
special_bonds   lj/coul 0.0 0.0 0.5

# ========================================
# Settings
# ========================================
neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes

# ========================================
# Energy Minimization
# ========================================
minimize        1.0e-4 1.0e-6 1000 10000

# ========================================
# Equilibration (NPT to get correct density)
# ========================================
velocity        all create 300.0 482634
fix             1 all npt temp 300.0 300.0 100.0 iso 1.0 1.0 1000.0

timestep        1.0
thermo          1000
thermo_style    custom step temp pe ke etotal press density vol

run             50000    # 50 ps equilibration

# ========================================
# Production (NVT)
# ========================================
unfix           1
fix             2 all nvt temp 300.0 300.0 100.0

reset_timestep  0
dump            1 all custom 1000 methanol.lammpstrj id mol type q x y z
dump_modify     1 sort id

# Compute RDF (O-O, shows H-bonding)
compute         rdf_oo all rdf 100 2 2
fix             3 all ave/time 100 10 1000 c_rdf_oo[*] file methanol_rdf.dat mode vector

run             100000   # 100 ps production

write_data      methanol_final.data
```

:::{important} Common Error: "All bond coeffs are not set"
If you get this error during minimization, it means the number of bond/angle/dihedral types declared in the data file header doesn't match the coefficients defined in the LAMMPS input.

**Fix:** Make sure the data file builder declares:
- `3 bond types` (not 4)
- `3 angle types` (not 5)
- `1 dihedral types` (not 2)

Re-run the `build_methanol.ipynb` notebook to regenerate `methanol.data` with correct type counts.
:::

### 1.5 Analysis: Hydrogen Bonding

**Using pre-run results:**

```bash
ls /anvil/projects/x-chm250117/class_examples/tutorial12_completed/part1_opls_methanol/methanol_rdf.dat
ls /anvil/projects/x-chm250117/class_examples/tutorial12_completed/part1_opls_methanol/methanol.log
```

Create `analysis.ipynb`:

**Notebook Cell 1: RDF Analysis**

```python
import numpy as np
import matplotlib.pyplot as plt

# Read the RDF file - it contains multiple output blocks
# We want the LAST block (final time-averaged data)
with open('/anvil/projects/x-chm250117/class_examples/tutorial12_completed/part1_opls_methanol/methanol_rdf.dat', 'r') as f:
    lines = f.readlines()

# Find all timestep markers (lines with just two numbers)
data_blocks = []
current_block = []

for line in lines:
    line = line.strip()
    if line.startswith('#'):
        continue
    
    parts = line.split()
    if len(parts) == 2:  # Timestep marker: "timestep num_bins"
        if current_block:
            data_blocks.append(current_block)
        current_block = []
    elif len(parts) == 4:  # Data line: "bin r g(r) coord"
        current_block.append([float(x) for x in parts])

# Add last block
if current_block:
    data_blocks.append(current_block)

# Use the last (most recent) block
data = np.array(data_blocks[-1])
r = data[:, 1]      # Column 2: r (distance)
g_r = data[:, 2]    # Column 3: g(r)

plt.figure(figsize=(10, 6))
plt.plot(r, g_r, 'b-', linewidth=2.5, label='OPLS-AA Methanol')
plt.axhline(y=1.0, color='k', linestyle=':', alpha=0.5)
plt.axvline(x=2.8, color='red', linestyle='--', alpha=0.7, 
            label='H-bond distance (~2.8 Å)')

plt.xlabel('O-O Distance (Å)', fontsize=12)
plt.ylabel('g(r)', fontsize=12)
plt.title('Methanol O-O RDF: Hydrogen Bonding Structure', fontsize=13, fontweight='bold')
plt.xlim(0.5, 10)
plt.ylim(0, 4)
plt.legend(fontsize=11)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('methanol_rdf.png', dpi=150)
plt.show()

# Find first peak (H-bonding)
first_peak_idx = np.argmax(g_r[:60])
print(f"\nFirst peak at r = {r[first_peak_idx]:.2f} Å")
print(f"Peak height: {g_r[first_peak_idx]:.2f}")
print(f"\nNumber of data blocks found: {len(data_blocks)}")
print(f"Using final block with {len(data)} bins")
print("\n→ First peak represents H-bonded neighbors")
print("→ OPLS-AA captures methanol H-bonding structure")
```

:::{tip} Understanding LAMMPS RDF Output
The `methanol_rdf.dat` file format:
```
# Time-averaged data for fix 3
# TimeStep Number-of-rows
# Row c_rdf_oo[1] c_rdf_oo[2] c_rdf_oo[3]
1000 100
1 2.05 0.00 0
2 2.15 0.12 45
...
```

- **Line 1-3:** Header comments (skip with `skiprows=4`)
- **Line 4:** Timestep and number of bins
- **Data columns:**
  - Column 1: Bin number
  - Column 2: Distance r (Å)
  - Column 3: g(r) - radial distribution function
  - Column 4: Coordination number (cumulative integral)

The `ave/time` fix averages over multiple frames (100 samples × 10 frequency = last 1000 steps).
:::

**Notebook Cell 2: Density Check**

```python
import re

# Set path to pre-run simulation
sim_dir = "/anvil/projects/x-chm250117/class_examples/tutorial12_completed/part1_opls_methanol"

# Extract density from log file
with open(f'{sim_dir}/methanol.log', 'r') as f:
    log_text = f.read()

# Find density values (last 10000 steps)
density_pattern = r'^\s*\d+\s+[\d.]+\s+[\d.e+-]+\s+[\d.e+-]+\s+[\d.e+-]+\s+[\d.e+-]+\s+([\d.]+)'
densities = []

for line in log_text.split('\n'):
    match = re.match(density_pattern, line)
    if match:
        densities.append(float(match.group(1)))

if len(densities) > 100:
    # Take last 100 values (production run)
    production_density = densities[-100:]
    mean_rho = np.mean(production_density)
    std_rho = np.std(production_density)
    
    print(f"Density Analysis:")
    print(f"  Mean: {mean_rho:.4f} g/cm³")
    print(f"  Std:  {std_rho:.4f} g/cm³")
    print(f"  Experimental (298K): 0.791 g/cm³")
    print(f"  Error: {100*(mean_rho - 0.791)/0.791:.1f}%")
    print("\n→ OPLS-AA designed for liquid phase properties")
    print("→ Should match experimental density within ~5%")
else:
    print("Not enough data points. Check log file.")
```

:::{note} Force Field Validation
The key test of a force field is whether it reproduces experimental properties:
- **Density**: Should be ~0.79 g/cm³ at 300K
- **RDF first peak**: ~2.8 Å (H-bonding)
- **Heat of vaporization**: ~37.5 kJ/mol (requires longer runs)

If your results are far off, check:
1. Parameter values were typed correctly
2. Equilibration was long enough
3. System size is sufficient (125 molecules minimum)
:::

---

## Part 2: EAM for Metals (Copper Nanoparticle)

### 2.1 Why EAM for Metals

Pair potentials (LJ) fail for metals because:
- Cohesive energy depends on coordination number
- Surface atoms bind differently than bulk
- Elastic constants cannot be fit correctly

**EAM Solution:** Bonding energy depends on local electron density.

### 2.2 System: Nanoparticle Surface Melting

We'll create a **spherical copper nanoparticle** (~1000 atoms, ~4 nm diameter) with **free surfaces**.

**Key difference from bulk:**
- **Shrink-wrapped boundaries** (`boundary s s s`) = isolated finite particle
- **Surface atoms** have fewer neighbors → weaker binding
- **Melts 100-200 K below bulk** Tm (surface melting)
- **Expected melting range:** ~1150-1250 K (vs. bulk 1358 K)

### 2.3 LAMMPS Input File

The LAMMPS input file `copper_melt.in`:

```lammps
# Copper Nanoparticle Melting - EAM Potential
# Tutorial 12 - Part 2
# Shrink-wrapped boundaries = isolated nanoparticle with free surfaces

units           metal
atom_style      atomic
boundary        s s s          # Shrink-wrapped = non-periodic (finite particle)

# ========================================
# Create Spherical Copper Nanoparticle
# ========================================
lattice         fcc 3.615
region          box block 0 10 0 10 0 10
create_box      1 box

# Create sphere in center of box (radius in lattice units)
region          sphere sphere 5 5 5 3.5 units lattice
create_atoms    1 region sphere

mass            1 63.546

print           "Created nanoparticle with $(count(all)) atoms"

# ========================================
# EAM Potential
# ========================================
pair_style      eam
pair_coeff      1 1 /apps/spack/anvil/apps/lammps/20210310-gcc-11.2.0-jzfe7x3/share/lammps/potentials/Cu_u3.eam

# ========================================
# Settings
# ========================================
neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes

# ========================================
# Equilibration at 300K (solid nanoparticle)
# ========================================
velocity        all create 300.0 857321
fix             1 all nvt temp 300.0 300.0 0.1

timestep        0.001
thermo          500
thermo_style    custom step temp pe ke etotal

run             5000           # 5 ps equilibration

# ========================================
# Heating: 300K → 1400K (nanoparticle melts at lower T)
# ========================================
# Nanoparticles melt 100-200K below bulk due to surface effects

unfix           1
fix             2 all nvt temp 300.0 1400.0 0.1

# Compute coordination number
compute         coord all coord/atom cutoff 3.2

# Compute MSD
compute         msd all msd

# Output
dump            1 all custom 500 copper_heat.lammpstrj id type x y z c_coord
dump_modify     1 sort id

fix             3 all ave/time 50 1 50 c_msd[4] file msd.dat

thermo_style    custom step temp pe ke etotal c_msd[4]

run             100000         # 100 ps heating

write_data      copper_final.data
```

### 2.4 Locating the EAM Potential File

:::{important} EAM File Location
LAMMPS comes with many EAM potential files. On Anvil, they're in:
```
/apps/spack/anvil/apps/lammps/20210310-gcc-11.2.0-jzfe7x3/share/lammps/potentials/
```

The tutorial uses the full path directly in the `pair_coeff` command, so you don't need to copy files.

Common EAM files:
- `Cu_u3.eam` - Copper (Mishin 2001)
- `Al_u3.eam` - Aluminum  
- `Ni_u3.eam` - Nickel
- `FeNiCr.eam.alloy` - Steel alloys (use `pair_style eam/alloy` for multi-element)
:::

### 2.5 Analysis: Detecting Melting

**Using pre-run results:**

```bash
ls /anvil/projects/x-chm250117/class_examples/tutorial12_completed/part2_eam_copper/copper_melt.log
ls /anvil/projects/x-chm250117/class_examples/tutorial12_completed/part2_eam_copper/copper_heat.lammpstrj
```

Create `analysis.ipynb`:

**Notebook Cell 1: Temperature and Energy vs. Time**

```python
import numpy as np
import matplotlib.pyplot as plt

# Set path to pre-run simulation
sim_dir = "/anvil/projects/x-chm250117/class_examples/tutorial12_completed/part2_eam_copper"

# Parse LAMMPS log file handling variable column counts
def parse_lammps_log(filename):
    """Extract thermo data from LAMMPS log file"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    data_blocks = []
    current_block = []
    in_data = False
    current_ncols = None
    
    for line in lines:
        stripped = line.strip()
        
        # Start of thermo output
        if stripped.startswith('Step'):
            if current_block:  # Save previous block
                data_blocks.append(np.array(current_block))
                current_block = []
            in_data = True
            current_ncols = len(stripped.split())
            continue
        
        # End of thermo output
        if in_data and (stripped.startswith('Loop') or 
                        stripped == '' or
                        stripped.startswith('WARNING') or
                        stripped.startswith('Performance')):
            in_data = False
            if current_block:
                data_blocks.append(np.array(current_block))
                current_block = []
            continue
        
        # Collect data lines
        if in_data:
            try:
                values = [float(x) for x in stripped.split()]
                if len(values) == current_ncols:
                    current_block.append(values)
            except (ValueError, IndexError):
                continue
    
    # Save last block
    if current_block:
        data_blocks.append(np.array(current_block))
    
    return data_blocks

# Load the data
blocks = parse_lammps_log(f'{sim_dir}/copper_melt.log')

print(f"Found {len(blocks)} thermo output blocks")
for i, block in enumerate(blocks):
    print(f"  Block {i+1}: {len(block)} rows × {block.shape[1]} columns")

# Block 1 = equilibration (5 cols: Step Temp PE KE TotEng)
# Block 2 = heating (6 cols: Step Temp PE KE TotEng MSD)

if len(blocks) < 2:
    print("Error: Expected 2 thermo blocks (equilibration + heating)")
else:
    heating_data = blocks[1]  # Second block has MSD
    
    temps = heating_data[:, 1]  # Temperature
    pe = heating_data[:, 2]     # Potential Energy  
    msd = heating_data[:, 5]    # MSD (last column)
    
    # Create plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: Potential Energy vs. Temperature
    ax1.plot(temps, pe, 'b-', linewidth=2)
    ax1.axvline(x=1358, color='red', linestyle='--', linewidth=1.5, 
                alpha=0.7, label='Bulk Cu Tm = 1358 K')
    ax1.set_xlabel('Temperature (K)', fontsize=12)
    ax1.set_ylabel('Potential Energy (eV)', fontsize=12)
    ax1.set_title('PE vs. T (Plateau = latent heat)', fontsize=13, fontweight='bold')
    ax1.legend()
    ax1.grid(alpha=0.3)
    
    # Plot 2: MSD vs. Temperature
    ax2.plot(temps, msd, 'r-', linewidth=2)
    ax2.axvline(x=1358, color='blue', linestyle='--', linewidth=1.5, 
                alpha=0.7, label='Expected Tm')
    ax2.set_xlabel('Temperature (K)', fontsize=12)
    ax2.set_ylabel('MSD (Å²)', fontsize=12)
    ax2.set_title('MSD vs. T (Jump = liquid diffusion)', fontsize=13, fontweight='bold')
    ax2.set_yscale('log')
    ax2.legend()
    ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('copper_melting.png', dpi=150)
    plt.show()
    
    # Estimate melting point from MSD jump
    msd_threshold = 5.0  # Lower threshold for nanoparticle
    melting_indices = np.where(msd > msd_threshold)[0]
    if len(melting_indices) > 0:
        approx_tm = temps[melting_indices[0]]
        print(f"\nApproximate melting point: {approx_tm:.0f} K")
        print(f"Expected (bulk Cu): 1358 K")
        print(f"Expected (4 nm nanoparticle): ~1150-1250 K")
        print(f"Depression: {1358 - approx_tm:.0f} K below bulk")
        print("\nNote: Surface melting occurs at lower T than bulk")
    else:
        print(f"\nNo clear melting detected (MSD max = {msd.max():.1f} Å²)")
        print("Try heating to higher temperature")
```

**Notebook Cell 2: Coordination Number (Structure Change)**

```python
print("="*60)
print("RECOMMENDED: Visualize in OVITO")
print("="*60)
print(f"\nTrajectory location:")
print(f"  /anvil/projects/x-chm250117/class_examples/tutorial12_completed/part2_eam_copper/copper_heat.lammpstrj")
print("\n1. Open: copper_heat.lammpstrj in OVITO")
print("2. Add modification → Structure identification → Coordination analysis")
print("3. Cutoff: 3.2 Å (first neighbor shell)")
print("4. Plot coordination vs. frame number")
print("\nExpected:")
print("  Solid core: ~12 neighbors (bulk-like)")
print("  Surface atoms: ~6-9 neighbors (undercoordinated)")
print("  Liquid: ~10-11 neighbors globally")
print("  Surface melts first → core melts later")

print("\nManual coordination check:")
print("  Solid FCC Cu: coordination ~12 (bulk), ~9 (surface)")
print("  Liquid Cu: coordination ~11 (more disordered)")
```

:::{tip} Observing Melting in VMD
To visualize the melting transition:
1. Load `copper_heat.lammpstrj` in VMD
2. Color atoms by coordination number (c_coord)
3. Play animation through heating
4. You'll see: solid (ordered) → liquid (disordered)
:::

---

## Part 3: ReaxFF Reactive Chemistry (Methane Combustion)

### 3.1 The Goal: Simulating Chemical Reactions

ReaxFF allows bonds to break and form dynamically. We'll simulate:
```
CH₄ + 2O₂ → CO₂ + 2H₂O
```

At high temperature (4000 K), methane combusts without pre-defined reaction pathways.

### 3.2 Building the System

The data file `combustion.data` is already provided in the shared directory:

```bash
ls /anvil/projects/x-chm250117/class_examples/tutorial12_completed/part3_reaxff_combustion/combustion.data
```

:::{dropdown} Python Script to Build Combustion System (for reference)
```python
import numpy as np

# System parameters
n_methane = 10
n_oxygen = 20      # 2:1 ratio for stoichiometric combustion
box_size = 15.0    # Smaller box for higher density, faster reactions

print("Building CH4 + O2 combustion system...")
print(f"  {n_methane} CH4 molecules")
print(f"  {n_oxygen} O2 molecules")
print(f"  Total atoms: {n_methane*5 + n_oxygen*2}")

# Methane geometry
ch_bond = 1.09
theta_tet = np.radians(109.47)

# CH4: C at center, H in tetrahedral positions
def make_methane(center):
    c = center
    h1 = c + ch_bond * np.array([1, 0, 0])
    h2 = c + ch_bond * np.array([-1/3, 2*np.sqrt(2)/3, 0])
    h3 = c + ch_bond * np.array([-1/3, -np.sqrt(2)/3, np.sqrt(2/3)])
    h4 = c + ch_bond * np.array([-1/3, -np.sqrt(2)/3, -np.sqrt(2/3)])
    return [c, h1, h2, h3, h4]

# O2: Simple dumbbell
def make_o2(center):
    o1 = center + np.array([0.6, 0, 0])
    o2 = center - np.array([0.6, 0, 0])
    return [o1, o2]

with open('combustion.data', 'w') as f:
    f.write('ReaxFF Combustion: CH4 + O2\n\n')
    
    n_atoms = n_methane * 5 + n_oxygen * 2
    f.write(f'{n_atoms} atoms\n')
    f.write('0 bonds\n')      # ReaxFF doesn't use topology
    f.write('0 angles\n')
    f.write('0 dihedrals\n\n')
    
    f.write('3 atom types  # 1=C, 2=H, 3=O\n\n')
    
    f.write(f'0.0 {box_size} xlo xhi\n')
    f.write(f'0.0 {box_size} ylo yhi\n')
    f.write(f'0.0 {box_size} zlo zhi\n\n')
    
    f.write('Masses\n\n')
    f.write('1 12.011  # C\n')
    f.write('2 1.008   # H\n')
    f.write('3 15.999  # O\n\n')
    
    f.write('Atoms  # charge\n\n')
    
    atom_id = 1
    
    # Place CH4 molecules randomly
    for i in range(n_methane):
        center = np.random.rand(3) * box_size
        atoms = make_methane(center)
        types = [1, 2, 2, 2, 2]
        
        for pos, atype in zip(atoms, types):
            f.write(f'{atom_id} {atype} 0.0 {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n')
            atom_id += 1
    
    # Place O2 molecules randomly
    for i in range(n_oxygen):
        center = np.random.rand(3) * box_size
        atoms = make_o2(center)
        
        for pos in atoms:
            f.write(f'{atom_id} 3 0.0 {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n')
            atom_id += 1

print(f"\n✓ Created combustion.data")
print(f"✓ Total atoms: {n_atoms}")
print("✓ Ready for ReaxFF simulation!")

import os
if os.path.exists('combustion.data'):
    size = os.path.getsize('combustion.data') / 1024
    print(f"✓ File size: {size:.1f} KB")
```
:::

### 3.3 LAMMPS Input with ReaxFF

The LAMMPS input file `combustion.in`:

```lammps
# Methane Combustion with ReaxFF
# Tutorial 12 - Part 3

units           real
atom_style      charge
boundary        p p p

read_data       combustion.data

# ========================================
# ReaxFF Potential
# ========================================
# CHO force field (hydrocarbons + oxygen)
pair_style      reax/c NULL safezone 3.0 mincap 150
pair_coeff      * * /apps/spack/anvil/apps/lammps/20210310-gcc-11.2.0-jzfe7x3/share/lammps/potentials/ffield.reax.cho C H O

# ReaxFF requires charge equilibration
fix             qeq all qeq/reax 1 0.0 10.0 1.0e-6 reax/c

# ========================================
# Settings
# ========================================
neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes

# ========================================
# Energy Minimization
# ========================================
minimize        1.0e-4 1.0e-6 1000 10000

# ========================================
# Heat to Combustion Temperature
# ========================================
velocity        all create 300.0 458293
fix             1 all nvt temp 300.0 4000.0 100.0

timestep        0.1            # 0.1 fs (ReaxFF needs small timestep)
thermo          100
thermo_style    custom step temp pe ke etotal press

run             10000          # Heat up

# ========================================
# Combustion Simulation at 4000K
# ========================================
unfix           1
fix             2 all nvt temp 4000.0 4000.0 100.0

reset_timestep  0

# Track species (bond order analysis)
fix             3 all reax/c/species 100 1 100 species.out element C H O

dump            1 all custom 100 combustion.lammpstrj id type q x y z
dump_modify     1 sort id

run             50000          # 50 ps at 4000K (reactions occur)

write_data      combustion_final.data
```

### 3.4 ReaxFF Parameter File Location

:::{important} ReaxFF Force Field Files
ReaxFF parameters are system-specific. For CHO (hydrocarbons + oxygen), use:

**File:** `ffield.reax.cho`

**Full path on Anvil:**
```bash
/apps/spack/anvil/apps/lammps/20210310-gcc-11.2.0-jzfe7x3/share/lammps/potentials/ffield.reax.cho
```

**You don't need to copy this file** - the LAMMPS input already references the full path in the `pair_coeff` line.

**Other available ReaxFF files:**
- `ffield.reax.Fe_O_C_H` - Steel oxidation/corrosion
- `ffield.reax.Si_O_H` - Silica/glass
- `ffield.reax.CHON` - Explosives (RDX, HMX)

See full list:
```bash
ls /apps/spack/anvil/apps/lammps/20210310-gcc-11.2.0-jzfe7x3/share/lammps/potentials/ffield.reax.*
```
:::

### 3.5 Analysis: Visualizing Reactions in OVITO

**Using pre-run results:**

```bash
ls /anvil/projects/x-chm250117/class_examples/tutorial12_completed/part3_reaxff_combustion/combustion.lammpstrj
ls /anvil/projects/x-chm250117/class_examples/tutorial12_completed/part3_reaxff_combustion/species.out
```

Use OVITO to visually observe the combustion reactions.

**Open OVITO on Anvil:**
1. Go to Open OnDemand → Interactive Apps → OVITO
2. Launch OVITO session

**Load the trajectory:**
```
File → Load File → combustion.lammpstrj
```

**Recommended OVITO setup for visualizing reactions:**

**Step 1: Color atoms by type**
```
Add modification → Coloring → Color coding
- Property: Particle Type
- C (type 1): Gray or Black
- H (type 2): White
- O (type 3): Red
```

**Step 2: Adjust atom sizes**
```
Particles → Display settings:
- Radius: 0.4 Å (or adjust to taste)
- This makes atoms visible but allows seeing through structure
```

**Step 3: Create bonds to visualize molecules**
```
Add modification → Topology → Create bonds
- Cutoff radius: 1.8 Å
  (Captures C-C, C-H, C-O, O-H, O-O bonds)
- Show bonds: Enabled
- Bond width: 0.15 Å
```

**Step 4: Play the animation**
```
- Use timeline slider at bottom
- Or click play button
- Adjust speed with FPS control
```

**What to observe:**

**Initial state (frame 0):**
- ✓ Organized CH₄ molecules (1 C + 4 H atoms bonded)
- ✓ O₂ molecules (2 O atoms bonded)
- ✓ Clear separation between molecules

**During reaction (middle frames):**
- ✓ Bonds breaking (CH₄ dissociates, O₂ dissociates)
- ✓ Atoms mixing rapidly
- ✓ New bonds forming (C-O, O-H)
- ✓ Transient species (CO, H₂O forming)

**Final state (last frame):**
- ✓ CO₂ molecules (1 C + 2 O linear)
- ✓ H₂O molecules (1 O + 2 H bent)
- ✓ Possibly some CO, H₂ if incomplete combustion
- ✓ More organized than middle frames

:::{dropdown} Advanced OVITO Visualization Tips
**Color by potential energy per atom:**
```
Add modification → Coloring → Color coding
- Property: Potential Energy
- Gradient: Blue (low) → Red (high)
- Shows which atoms are in high-energy configurations
```

**Track coordination number:**
```
Add modification → Structure identification → Coordination analysis
- Cutoff: 1.8 Å
- Shows how many neighbors each atom has
- C should have ~3-4 (in CO₂, has 2 O neighbors)
- O in H₂O should have ~1-2 (1 H neighbor)
```

**Create snapshots:**
```
File → Export File → Image
- Save frames at: initial, middle, final
- Compare to see reaction progress
```

**Export bond statistics:**
```
Add modification → Analysis → Bond analysis
- File → Export File → Table
- Exports bond lengths and types over time
```
:::

:::{tip} Understanding What You See
**CH₄ → CO₂ + H₂O transformation:**

1. **CH₄ breakdown:**
   - 4 C-H bonds break
   - C atom freed
   - H atoms freed

2. **O₂ activation:**
   - O=O bond breaks
   - Reactive O atoms freed

3. **Product formation:**
   - C + 2O → O=C=O (linear CO₂)
   - 2H + O → H-O-H (bent water)

**What distinguishes CO₂ from other species:**
- **Linear geometry** (O-C-O angle = 180°)
- **Two C=O double bonds** (shorter than C-O single)
- **Symmetric** structure

**What distinguishes H₂O:**
- **Bent geometry** (H-O-H angle ≈ 104°)
- **Two O-H bonds**
- **Persistent structure** once formed
:::

**Alternative: Quick Python check for atom conservation**

If you want to verify atoms are conserved:

**Create `analysis.ipynb`:**

```python
import numpy as np

# Set path to pre-run simulation
sim_dir = "/anvil/projects/x-chm250117/class_examples/tutorial12_completed/part3_reaxff_combustion"

# Count atoms in first and last frame
def count_atoms_by_type(trajfile):
    """Count C, H, O atoms in first and last frames"""
    
    with open(trajfile, 'r') as f:
        lines = f.readlines()
    
    frames = []
    i = 0
    while i < len(lines):
        if 'ITEM: TIMESTEP' in lines[i]:
            frame_atoms = {'C': 0, 'H': 0, 'O': 0}
            
            # Skip to atoms section
            while i < len(lines) and 'ITEM: ATOMS' not in lines[i]:
                i += 1
            i += 1
            
            # Count atoms
            while i < len(lines) and not lines[i].startswith('ITEM:'):
                parts = lines[i].split()
                atom_type = int(parts[1])
                if atom_type == 1:  # C
                    frame_atoms['C'] += 1
                elif atom_type == 2:  # H
                    frame_atoms['H'] += 1
                elif atom_type == 3:  # O
                    frame_atoms['O'] += 1
                i += 1
            
            frames.append(frame_atoms)
        else:
            i += 1
    
    return frames[0], frames[-1]

first, last = count_atoms_by_type(f'{sim_dir}/combustion.lammpstrj')

print("Atom Conservation Check:")
print("="*40)
print(f"          C    H    O")
print(f"Initial:  {first['C']:2d}  {first['H']:2d}  {first['O']:2d}")
print(f"Final:    {last['C']:2d}  {last['H']:2d}  {last['O']:2d}")
print("="*40)

if first == last:
    print("✓ Atoms conserved (as expected)")
else:
    print("✗ ERROR: Atoms not conserved!")

# Calculate expected products from stoichiometry
initial_ch4 = 10
initial_o2 = 20

print(f"\nExpected reaction:")
print(f"  {initial_ch4} CH₄ + {2*initial_ch4} O₂ → {initial_ch4} CO₂ + {2*initial_ch4} H₂O")
print(f"\nAtom balance:")
print(f"  C: {initial_ch4} (in {initial_ch4} CO₂)")
print(f"  H: {4*initial_ch4} (in {2*initial_ch4} H₂O)")  
print(f"  O: {4*initial_ch4} (in {initial_ch4} CO₂ + {2*initial_ch4} H₂O)")
```

---

## Summary and Comparison

### What We Learned

| Potential Type | System | Key Feature | What It Predicts |
|----------------|--------|-------------|------------------|
| **Force Field (OPLS)** | Methanol | Validated parameters | Density, structure, H-bonding |
| **EAM** | Copper | Coordination-dependent | Melting point, surface energy |
| **ReaxFF** | Combustion | Dynamic bonds | Reaction products, mechanisms |

### Computational Cost Reality

From your simulations, compare timesteps/second:

| System | Timestep | Atoms | Performance | Notes |
|--------|----------|-------|-------------|-------|
| Methanol (OPLS) | 1 fs | 750 | Medium | Coulomb (PPM) dominates |
| Copper (EAM) | 1 fs | 2048 | Fast | No charges, local calculation |
| Combustion (ReaxFF) | 0.25 fs | 90 | **Very Slow** | Bond order + QEq expensive |

**Lesson:** ReaxFF is 10-100× slower than classical force fields. Use it only when chemistry is essential.

### Validation Checklist

Before trusting any advanced potential:

**For Force Fields:**
- ✓ Density within 5% of experiment
- ✓ RDF peaks match experimental/neutron data
- ✓ Used correct parameter set (OPLS vs AMBER vs CHARMM)

**For EAM:**
- ✓ Lattice constant matches experiment
- ✓ Cohesive energy reasonable
- ✓ Melting point within 100-200 K of experiment

**For ReaxFF:**
- ✓ Used correct force field (CHO, not SiO or FeCr)
- ✓ Products match expected chemistry
- ✓ Energy release reasonable (exothermic reactions)

---

## Troubleshooting Guide

### Common LAMMPS Errors

**"ERROR: Cannot open file ffield.reax.cho"**
- This error should NOT occur if you used the full path in pair_coeff
- Verify your LAMMPS input has the complete path:
  ```
  pair_coeff * * /apps/spack/anvil/apps/lammps/20210310-gcc-11.2.0-jzfe7x3/share/lammps/potentials/ffield.reax.cho C H O
  ```
- Check the file exists: `ls /apps/spack/anvil/apps/lammps/20210310-gcc-11.2.0-jzfe7x3/share/lammps/potentials/ffield.reax.cho`

**"ERROR: Illegal pair_coeff command"**
- Check atom type count matches data file
- Verify element order in pair_coeff matches data file

**"ERROR: Lost atoms" (ReaxFF)**
- Reduce timestep (try 0.1 fs instead of 0.25 fs)
- Increase neighbor list cutoff
- Check initial structure (no severe overlaps)

**"WARNING: Inconsistent image flags"**
- Box too small for cutoff
- Increase box size or reduce cutoff

### Performance Issues

If simulations are too slow:
- **ReaxFF**: Reduce system size (try 5 CH₄ + 10 O₂)
- **EAM**: Use smaller nanoparticle (4×4×4 unit cells)
- **OPLS**: Check PPM accuracy (try `pppm 1.0e-3` instead of `1.0e-4`)

---

## Exercise Questions

1. **OPLS Methanol:** How does the O-O RDF compare to water? Why is the first peak different?

2. **EAM Copper:** What is the simulated melting point? How does it compare to bulk copper (1358 K)? Why might it be different?

3. **ReaxFF:** Does the combustion go to completion? If not, what intermediates remain? What would happen at lower temperature (1500 K)?

4. **Force Field Comparison:** Run the same system (e.g., methanol) with OPLS and AMBER. Do densities differ? Why?

---

## Next Steps

You now have hands-on experience with the three major classes of MD potentials. In real research:

- Use **classical force fields** (OPLS/AMBER/CHARMM) for biomolecules, organic liquids
- Use **EAM** for metals, alloys, nanoparticles
- Use **ReaxFF** when bonds must break/form

**For your final projects**, choose the simplest potential that captures your physics!