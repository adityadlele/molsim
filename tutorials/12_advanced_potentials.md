# Tutorial 12: Advanced MD Potentials in LAMMPS

:::{warning} Tutorial Status
**These examples are not fully tested.** This note will be removed after validation. If you encounter errors, please report them to the instructor.
:::

**Objective:** Learn to use professional force field parameter files, simulate metallic systems with EAM, and run reactive chemistry simulations with ReaxFF.

---

## File Organization

Create a well-organized directory structure on Anvil:

```bash
cd $SCRATCH
mkdir tutorial12
cd tutorial12

# Create subdirectories for each part
mkdir part1_opls_methanol
mkdir part2_eam_copper
mkdir part3_reaxff_combustion
```

**Recommended structure:**
```
tutorial12/
├── part1_opls_methanol/
│   ├── build_methanol.ipynb    # Structure builder
│   ├── methanol.in             # LAMMPS input
│   ├── submit.sh               # SLURM script
│   └── analysis.ipynb          # Analysis
├── part2_eam_copper/
│   ├── copper_melt.in          # LAMMPS input
│   ├── submit.sh
│   └── analysis.ipynb
└── part3_reaxff_combustion/
    ├── build_methane_o2.ipynb  # Structure builder
    ├── combustion.in           # LAMMPS input
    ├── submit.sh
    └── analysis.ipynb
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

**Navigate to directory:**
```bash
cd $SCRATCH/tutorial12/part1_opls_methanol
```

**Create notebook: `build_methanol.ipynb`**

**Notebook Cell 1: Build Methanol Liquid Box**

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
    f.write('4 bond types  # 1=C-O, 2=C-H, 3=O-H, 4=dummy\n')
    f.write('5 angle types\n')
    f.write('2 dihedral types\n\n')
    
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

### 1.4 LAMMPS Input with OPLS Parameters

Create `methanol.in`:

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

### 1.5 Analysis: Hydrogen Bonding

Create `analysis.ipynb` after the simulation:

**Notebook Cell 1: RDF Analysis**

```python
import numpy as np
import matplotlib.pyplot as plt

# Load O-O RDF
data = np.loadtxt('methanol_rdf.dat')
r = data[:, 1]
g_r = data[:, 2]

plt.figure(figsize=(10, 6))
plt.plot(r, g_r, 'b-', linewidth=2.5, label='OPLS-AA Methanol')
plt.axhline(y=1.0, color='k', linestyle=':', alpha=0.5)
plt.axvline(x=2.8, color='red', linestyle='--', alpha=0.7, 
            label='H-bond distance (~2.8 Å)')

plt.xlabel('O-O Distance (Å)', fontsize=12)
plt.ylabel('g(r)', fontsize=12)
plt.title('Methanol O-O RDF: Hydrogen Bonding Structure', fontsize=13, fontweight='bold')
plt.xlim(2, 10)
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
print("\n→ First peak represents H-bonded neighbors")
print("→ OPLS-AA captures methanol H-bonding structure")
```

**Notebook Cell 2: Density Check**

```python
# Extract density from log file
import re

with open('log.lammps', 'r') as f:
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

### 2.2 System: Copper Melting

We'll simulate a copper nanoparticle and observe melting—a phenomenon that emerges naturally from EAM but not from pair potentials.

### 2.3 LAMMPS Input File

Navigate to directory:
```bash
cd $SCRATCH/tutorial12/part2_eam_copper
```

Create `copper_melt.in`:

```lammps
# Copper Nanoparticle Melting - EAM Potential
# Tutorial 12 - Part 2

units           metal          # eV, Angstroms, ps
atom_style      atomic
boundary        p p p

# ========================================
# Create FCC Copper Crystal
# ========================================
lattice         fcc 3.615      # Copper lattice constant (Angstrom)
region          box block 0 8 0 8 0 8
create_box      1 box
create_atoms    1 box

mass            1 63.546       # Copper atomic mass

# ========================================
# EAM Potential
# ========================================
# Uses pre-tabulated EAM file from LAMMPS distribution
pair_style      eam/alloy
pair_coeff      * * Cu_u3.eam Cu

# NOTE: Cu_u3.eam is a standard EAM potential for Cu
# Located in LAMMPS potentials directory
# On Anvil: /apps/spack/bell/apps/lammps/20210310-x2lj7m2/share/lammps/potentials/

# ========================================
# Settings
# ========================================
neighbor        0.3 bin
neigh_modify    delay 0 every 1 check yes

# ========================================
# Equilibration at 300K (solid)
# ========================================
velocity        all create 300.0 857321
fix             1 all nvt temp 300.0 300.0 0.1

timestep        0.001          # 1 fs
thermo          1000
thermo_style    custom step temp pe ke etotal press

run             10000          # 10 ps equilibration

# ========================================
# Heating: 300K → 1500K (through melting point)
# ========================================
# Copper melting point: ~1358 K

unfix           1
fix             2 all nvt temp 300.0 1500.0 0.1

# Compute coordination number (detect melting)
compute         coord all coord/atom cutoff 3.2  # First neighbor shell

# Compute MSD (detect diffusion)
compute         msd all msd

# Output
dump            1 all custom 1000 copper_heat.lammpstrj id type x y z c_coord
dump_modify     1 sort id

fix             3 all ave/time 100 1 100 c_msd[4] file msd.dat

thermo_style    custom step temp pe ke etotal press c_msd[4]

run             100000         # 100 ps heating

write_data      copper_final.data
```

### 2.4 Locating the EAM Potential File

:::{important} EAM File Location
LAMMPS comes with many EAM potential files. On Anvil, they're in:
```
/apps/spack/bell/apps/lammps/20210310-x2lj7m2/share/lammps/potentials/
```

To use them, either:
1. **Copy to your directory:**
   ```bash
   cp /apps/spack/bell/apps/lammps/20210310-x2lj7m2/share/lammps/potentials/Cu_u3.eam .
   ```

2. **Or use full path in LAMMPS input:**
   ```lammps
   pair_coeff * * /apps/spack/bell/apps/lammps/.../Cu_u3.eam Cu
   ```

Common EAM files:
- `Cu_u3.eam` - Copper (Mishin 2001)
- `Al_u3.eam` - Aluminum
- `Ni_u3.eam` - Nickel
- `FeNiCr.eam.alloy` - Steel alloys
:::

### 2.5 Analysis: Detecting Melting

Create `analysis.ipynb`:

**Notebook Cell 1: Temperature and Energy vs. Time**

```python
import numpy as np
import matplotlib.pyplot as plt

# Parse LAMMPS log file
with open('log.lammps', 'r') as f:
    lines = f.readlines()

# Find thermo output during heating run
temps = []
pe = []
msd = []

in_thermo = False
for line in lines:
    if 'Step' in line and 'Temp' in line:
        in_thermo = True
        continue
    if 'Loop time' in line:
        in_thermo = False
    
    if in_thermo and line.strip() and not line.startswith('#'):
        parts = line.split()
        if len(parts) >= 6:
            try:
                temps.append(float(parts[1]))
                pe.append(float(parts[2]))
                msd.append(float(parts[6]))
            except:
                pass

temps = np.array(temps)
pe = np.array(pe)
msd = np.array(msd)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: Potential energy vs. temperature
ax1.plot(temps, pe, 'b-', linewidth=2)
ax1.axvline(x=1358, color='red', linestyle='--', linewidth=2, 
            label='Experimental Tm = 1358 K')
ax1.set_xlabel('Temperature (K)', fontsize=12)
ax1.set_ylabel('Potential Energy (eV/atom)', fontsize=12)
ax1.set_title('Melting Transition', fontsize=13, fontweight='bold')
ax1.legend(fontsize=11)
ax1.grid(alpha=0.3)

# Plot 2: MSD vs. temperature (diffusion indicator)
ax2.plot(temps, msd, 'r-', linewidth=2)
ax2.axvline(x=1358, color='blue', linestyle='--', linewidth=2,
            label='Expected Tm')
ax2.set_xlabel('Temperature (K)', fontsize=12)
ax2.set_ylabel('Mean Squared Displacement (ų)', fontsize=12)
ax2.set_title('Diffusion Onset (Melting)', fontsize=13, fontweight='bold')
ax2.set_yscale('log')
ax2.legend(fontsize=11)
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('copper_melting.png', dpi=150)
plt.show()

print("\nMelting Point Detection:")
print("  Look for:")
print("  1. Slope change in PE vs. T (latent heat)")
print("  2. Sharp increase in MSD (liquid diffusion)")
print("  3. Should occur near 1358 K")
```

**Notebook Cell 2: Coordination Number (Structure Change)**

```python
# Load trajectory and analyze coordination
import MDAnalysis as mda

try:
    u = mda.Universe('copper_final.data', 'copper_heat.lammpstrj')
    
    # Extract coordination numbers over time
    coords = []
    for ts in u.trajectory[::10]:  # Every 10th frame
        # Get coordination from dump file (stored in vx column by LAMMPS)
        coord_vals = u.atoms.positions[:, 0]  # Placeholder - actual extraction varies
        coords.append(np.mean(coord_vals))
    
    print("Coordination analysis requires MDAnalysis parsing")
    print("Alternatively, parse dump file directly:")
    
except ImportError:
    print("MDAnalysis not installed. Parse dump file manually:")

# Alternative: Direct dump file parsing
print("\nManual coordination check:")
print("  Solid FCC Cu: coordination ~12 (bulk), ~9 (surface)")
print("  Liquid Cu: coordination ~11 (more disordered)")
print("  Use: grep 'c_coord' copper_heat.lammpstrj")
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

At high temperature (2500 K), methane combusts without pre-defined reaction pathways.

### 3.2 Building the System

Navigate to directory:
```bash
cd $SCRATCH/tutorial12/part3_reaxff_combustion
```

Create `build_methane_o2.ipynb`:

**Notebook Cell 1: Build Combustion System**

```python
import numpy as np

# System parameters
n_methane = 10
n_oxygen = 20      # 2:1 ratio for stoichiometric combustion
box_size = 25.0

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

### 3.3 LAMMPS Input with ReaxFF

Create `combustion.in`:

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
pair_coeff      * * ffield.reax.cho C H O

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
fix             1 all nvt temp 300.0 2500.0 100.0

timestep        0.25           # 0.25 fs (ReaxFF needs small timestep)
thermo          100
thermo_style    custom step temp pe ke etotal press

run             10000          # Heat up

# ========================================
# Combustion Simulation at 2500K
# ========================================
unfix           1
fix             2 all nvt temp 2500.0 2500.0 100.0

reset_timestep  0

# Track species (bond order analysis)
fix             3 all reax/c/species 100 1 100 species.out element C H O

dump            1 all custom 100 combustion.lammpstrj id type q x y z
dump_modify     1 sort id

run             50000          # 50 ps at 2500K (reactions occur)

write_data      combustion_final.data
```

### 3.4 Obtaining the ReaxFF Parameter File

:::{important} ReaxFF Force Field Files
ReaxFF parameters are system-specific. For CHO (hydrocarbons + oxygen), use:

**File:** `ffield.reax.cho`

**Location on Anvil:**
```bash
/apps/spack/bell/apps/lammps/20210310-x2lj7m2/share/lammps/potentials/ffield.reax.cho
```

**Copy to your directory:**
```bash
cp /apps/spack/bell/apps/lammps/.../ffield.reax.cho .
```

**Other available ReaxFF files:**
- `ffield.reax.Fe_O_C_H` - Steel oxidation/corrosion
- `ffield.reax.Si_O_H` - Silica/glass
- `ffield.reax.CHON` - Explosives (RDX, HMX)

See full list: `ls /apps/spack/.../potentials/ffield.reax.*`
:::

### 3.5 Analysis: Detecting Reactions

Create `analysis.ipynb`:

**Notebook Cell 1: Species Evolution**

```python
import numpy as np
import matplotlib.pyplot as plt
import re

# Parse species output file
print("Parsing ReaxFF species data...")

with open('species.out', 'r') as f:
    lines = f.readlines()

# Extract timestep and species counts
data = {}
current_step = None

for line in lines:
    # Look for timestep markers
    if line.startswith('#'):
        match = re.search(r'timestep\s+(\d+)', line)
        if match:
            current_step = int(match.group(1))
            if current_step not in data:
                data[current_step] = {}
    else:
        # Parse species lines: "CH4  5"
        parts = line.split()
        if len(parts) >= 2 and current_step is not None:
            species = parts[0]
            count = int(parts[1])
            data[current_step][species] = count

# Convert to arrays
steps = sorted(data.keys())
species_of_interest = ['CH4', 'O2', 'CO2', 'H2O', 'CO', 'H2']

species_counts = {sp: [] for sp in species_of_interest}

for step in steps:
    for sp in species_of_interest:
        species_counts[sp].append(data[step].get(sp, 0))

# Plot
plt.figure(figsize=(12, 6))

colors = {'CH4': 'blue', 'O2': 'red', 'CO2': 'green', 
          'H2O': 'cyan', 'CO': 'orange', 'H2': 'purple'}

for sp in species_of_interest:
    if max(species_counts[sp]) > 0:
        plt.plot(steps, species_counts[sp], 'o-', linewidth=2, 
                markersize=4, label=sp, color=colors.get(sp, 'gray'))

plt.xlabel('Timestep', fontsize=12)
plt.ylabel('Number of Molecules', fontsize=12)
plt.title('ReaxFF Combustion: Species Evolution', fontsize=14, fontweight='bold')
plt.legend(fontsize=11, loc='best')
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('species_evolution.png', dpi=150)
plt.show()

print("\nExpected Reaction:")
print("  CH₄ + 2O₂ → CO₂ + 2H₂O")
print("\nObserve:")
print("  • CH4 and O2 decrease")
print("  • CO2 and H2O increase")
print("  • Possible intermediates: CO, H2, radicals")
```

**Notebook Cell 2: Final Composition**

```python
# Final state analysis
final_step = steps[-1]
final_composition = data[final_step]

print(f"\nFinal Composition (step {final_step}):")
print("="*40)

total_species = 0
for species, count in sorted(final_composition.items(), key=lambda x: -x[1]):
    if count > 0:
        print(f"  {species:<10} {count:>3}")
        total_species += count

print("="*40)
print(f"Total molecules: {total_species}")

# Check stoichiometry
initial_ch4 = species_counts['CH4'][0]
initial_o2 = species_counts['O2'][0]
final_co2 = species_counts['CO2'][-1]
final_h2o = species_counts['H2O'][-1]

print(f"\nReaction Progress:")
print(f"  Initial: {initial_ch4} CH₄ + {initial_o2} O₂")
print(f"  Final:   {final_co2} CO₂ + {final_h2o} H₂O")

# Expected products for complete combustion
expected_co2 = initial_ch4
expected_h2o = initial_ch4 * 2

print(f"  Expected: {expected_co2} CO₂ + {expected_h2o} H₂O")

conversion = 100 * (initial_ch4 - species_counts['CH4'][-1]) / initial_ch4
print(f"\nCH₄ Conversion: {conversion:.1f}%")

if conversion > 80:
    print("→ Good conversion! Most CH₄ reacted")
elif conversion > 50:
    print("→ Partial conversion. Try longer simulation or higher T")
else:
    print("→ Low conversion. Check temperature and mixing")
```

:::{note} Understanding ReaxFF Output
The `species.out` file tracks molecules identified by bond-order analysis:
- **Timestep headers**: `# Timestep XXXX`
- **Species lines**: `CH4 5` means 5 methane molecules detected
- Bond orders > 0.3 are considered "bonds"
- Species are identified by connectivity, not fixed topology

Common intermediates in combustion:
- **CO** - carbon monoxide (incomplete combustion)
- **H2** - hydrogen gas
- **OH, H, O** - radical species (very short-lived)
:::

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
- Solution: Copy force field file to working directory
- Check path: `ls ffield.reax.cho`

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