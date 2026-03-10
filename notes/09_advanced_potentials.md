# Advanced MD Potentials: Force Fields, Metals, and Reactive Chemistry

---

## 1. Introduction: Beyond Simple Molecules

In the previous lecture, we built up the classical force field piece by piece:
- LJ for van der Waals forces
- Harmonic bonds, angles, and dihedrals for molecular structure
- Coulomb interactions for electrostatics

These components form the foundation of **molecular mechanics force fields**—complete, self-consistent parameter sets designed to simulate specific classes of molecules.

However, there are two major classes of materials where classical force fields fail completely:

1. **Metals and Alloys**: Metallic bonding involves delocalized electrons in a "sea" that cannot be described by pair potentials or fixed bonds.

2. **Reactive Chemistry**: Bond breaking and forming (combustion, catalysis, polymerization) requires bonds to appear and disappear dynamically during the simulation.

This lecture covers three advanced topics:
- **Complete Force Fields** (AMBER, CHARMM, OPLS) - How professional parameter sets are organized
- **Many-Body Potentials** (EAM) - Metallic systems
- **Reactive Potentials** (ReaxFF) - Chemical reactions

---

## 2. Complete Biomolecular Force Fields

### 2.1 The Philosophy: Atom Typing

In classical force fields, not all carbon atoms are treated the same. A carbon in a benzene ring has different electronic structure (sp² hybridization) than a carbon in methane (sp³ hybridization). This leads to the concept of **atom types**.

**Atom Type** = Chemical element + Local bonding environment

**Example: Carbon atom types in AMBER:**
- `C` - sp² carbonyl carbon (C=O)
- `CA` - sp² aromatic carbon (benzene)
- `CT` - sp³ aliphatic carbon (saturated, like CH₃)
- `C*` - sp² carbon in 5-membered aromatic ring

Each atom type has its own set of parameters (mass, charge, LJ σ and ε).

### 2.2 The Three Major Force Fields for Biomolecules

#### AMBER (Assisted Model Building with Energy Refinement)

**History:** Developed by Peter Kollman's group at UCSF starting in the 1980s.

**Strengths:**
- Excellent for proteins and nucleic acids
- Well-validated for biomolecular simulations in water
- Many variants: ff14SB (proteins), ff19SB (latest), OL15 (RNA)

**Functional Form:**
$$U_{AMBER} = \sum_{bonds} k_b(r-r_0)^2 + \sum_{angles} k_\theta(\theta-\theta_0)^2$$
$$+ \sum_{dihedrals} \frac{V_n}{2}[1+\cos(n\phi-\gamma)] + \sum_{i<j}\left[\frac{A_{ij}}{r_{ij}^{12}} - \frac{B_{ij}}{r_{ij}^6} + \frac{q_iq_j}{4\pi\epsilon_0 r_{ij}}\right]$$

**Key Feature:** Uses **RESP** (Restrained Electrostatic Potential) charges—partial charges fitted to reproduce the electrostatic potential from quantum calculations.

#### CHARMM (Chemistry at Harvard Macromolecular Mechanics)

**History:** Developed by Martin Karplus's group at Harvard starting in the 1980s.

**Strengths:**
- Very comprehensive (proteins, lipids, carbohydrates, small molecules)
- CHARMM-GUI provides excellent web-based system builder
- Widely used for membrane proteins

**Functional Form:** Similar to AMBER but includes **Urey-Bradley terms** (1-3 distance restraints) and **CMAP** (backbone correction maps):

$$U_{CHARMM} = U_{bonds} + U_{angles} + U_{Urey-Bradley} + U_{dihedrals} + U_{impropers} + U_{CMAP} + U_{non-bonded}$$

**Urey-Bradley term:**
$$U_{UB} = k_{ub}(r_{13} - r_{ub,0})^2$$
This restrains the distance between atoms 1 and 3 in an angle (1-2-3), providing extra rigidity.

**CMAP (Correction Map):**
A 2D grid correction applied to protein backbone (φ, ψ) dihedral pairs to improve secondary structure prediction.

#### OPLS (Optimized Potentials for Liquid Simulations)

**History:** Developed by William Jorgensen's group at Yale.

**Strengths:**
- Specifically optimized for **liquid-phase** properties
- Excellent for small organic molecules, solvents
- OPLS-AA (all-atom) widely used in drug design

**Functional Form:** Similar to AMBER, but:
- Different combining rules for LJ parameters (geometric mean for both σ and ε)
- 1-4 scaling factors: 0.5 for both LJ and Coulomb

**Philosophy:** Parameters fitted to reproduce **condensed-phase** properties (liquid density, heat of vaporization) rather than gas-phase quantum calculations.

### 2.3 Interactive Example: Comparing Atom Types

Let's visualize how different atom types affect LJ parameters:

```python
import numpy as np
import matplotlib.pyplot as plt

# Example atom types from AMBER force field
atom_types = {
    'H': {'sigma': 0.6000, 'epsilon': 0.0157, 'description': 'Hydrogen (non-polar)'},
    'HC': {'sigma': 1.4870, 'epsilon': 0.0157, 'description': 'Hydrogen (aliphatic C-H)'},
    'HO': {'sigma': 0.0000, 'epsilon': 0.0000, 'description': 'Hydrogen (hydroxyl O-H)'},
    'CT': {'sigma': 3.3996, 'epsilon': 0.1094, 'description': 'Carbon (sp³ aliphatic)'},
    'CA': {'sigma': 3.3996, 'epsilon': 0.0860, 'description': 'Carbon (sp² aromatic)'},
    'C': {'sigma': 3.3996, 'epsilon': 0.0860, 'description': 'Carbon (sp² carbonyl)'},
    'N': {'sigma': 3.2500, 'epsilon': 0.1700, 'description': 'Nitrogen (amide)'},
    'O': {'sigma': 2.9599, 'epsilon': 0.2100, 'description': 'Oxygen (carbonyl)'},
    'OH': {'sigma': 3.0665, 'epsilon': 0.2104, 'description': 'Oxygen (hydroxyl)'},
}

# Create LJ potential curves for selected atom types
r = np.linspace(2.5, 8, 300)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# --- Plot 1: Carbon atom types ---
carbon_types = ['CT', 'CA', 'C']
colors_c = ['blue', 'green', 'purple']

for atype, color in zip(carbon_types, colors_c):
    sigma = atom_types[atype]['sigma']
    epsilon = atom_types[atype]['epsilon']
    U_LJ = 4 * epsilon * ((sigma/r)**12 - (sigma/r)**6)
    ax1.plot(r, U_LJ, color=color, linewidth=2.5, 
             label=f"{atype}: {atom_types[atype]['description']}")

ax1.axhline(y=0, color='k', linestyle='--', alpha=0.3)
ax1.set_xlabel('Distance r (Å)', fontsize=12)
ax1.set_ylabel('LJ Energy (kcal/mol)', fontsize=12)
ax1.set_title('Carbon Atom Types in AMBER', fontsize=13, fontweight='bold')
ax1.set_xlim(2.5, 8)
ax1.set_ylim(-0.15, 0.5)
ax1.legend(fontsize=10)
ax1.grid(alpha=0.3)

# --- Plot 2: Heteroatoms ---
hetero_types = ['N', 'O', 'OH']
colors_h = ['orange', 'red', 'darkred']

for atype, color in zip(hetero_types, colors_h):
    sigma = atom_types[atype]['sigma']
    epsilon = atom_types[atype]['epsilon']
    U_LJ = 4 * epsilon * ((sigma/r)**12 - (sigma/r)**6)
    ax2.plot(r, U_LJ, color=color, linewidth=2.5,
             label=f"{atype}: {atom_types[atype]['description']}")

ax2.axhline(y=0, color='k', linestyle='--', alpha=0.3)
ax2.set_xlabel('Distance r (Å)', fontsize=12)
ax2.set_ylabel('LJ Energy (kcal/mol)', fontsize=12)
ax2.set_title('Heteroatom Types in AMBER', fontsize=13, fontweight='bold')
ax2.set_xlim(2.5, 8)
ax2.set_ylim(-0.25, 0.5)
ax2.legend(fontsize=10)
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('atom_types_comparison.png', dpi=150, bbox_inches='tight')
plt.show()

print("\n" + "="*60)
print("ATOM TYPE PARAMETER COMPARISON:")
print("="*60)
print(f"{'Type':<6} {'σ (Å)':<8} {'ε (kcal/mol)':<15} Description")
print("-" * 60)
for atype, params in atom_types.items():
    print(f"{atype:<6} {params['sigma']:<8.4f} {params['epsilon']:<15.4f} {params['description']}")
print("="*60)
print("\nKey observations:")
print("  • Aromatic C (CA) has weaker ε than aliphatic C (CT)")
print("  • Oxygen types are smaller (σ) but stronger (ε)")
print("  • Hydrogen atoms vary widely - some have zero LJ!")
print("="*60)
```

### 2.4 Parameter File Structure

Force fields are distributed as **parameter files** that contain:

1. **Atom types and masses**
2. **Bond parameters** (k_b, r_0 for each bond type)
3. **Angle parameters** (k_θ, θ_0 for each angle type)
4. **Dihedral parameters** (V_n, n, γ for each dihedral type)
5. **Non-bonded parameters** (σ, ε for each atom type)
6. **Combining rules** (how to mix different atom types)

**Example snippet from AMBER parm10.dat:**
```
CT-CT  310.0    1.526       ! Aliphatic C-C bond
CT-HC  340.0    1.090       ! Aliphatic C-H bond
C -O   570.0    1.229       ! Carbonyl C=O bond

CT-CT-CT   40.0      109.50   ! C-C-C angle (tetrahedral)
HC-CT-HC   35.0      109.50   ! H-C-H angle

CT-CT-CT-CT  1  0.18  0.0 -3.0  1.0   ! Alkane torsion (3-fold)
```

**Critical:** Never modify these files unless you know exactly what you're doing. Force fields are internally consistent—changing one parameter breaks the balance.

### 2.5 Water Models: A Special Case

Water is so important in biomolecular simulations that it has its own family of models:

| Model | Sites | LJ Sites | Charges | Geometry | Comments |
|-------|-------|----------|---------|----------|----------|
| **TIP3P** | 3 | O only | 3 | Rigid | Fast, standard in AMBER/CHARMM |
| **TIP4P** | 4 | O only | 4 (virtual site) | Rigid | Better structure than TIP3P |
| **TIP5P** | 5 | O only | 5 (2 virtual) | Rigid | Best density, tetrahedral |
| **SPC/E** | 3 | O only | 3 | Rigid | Standard in GROMACS |
| **OPC** | 3 | O only | 3 | Rigid | Best all-around (2014) |

```python
# Compare water models
water_models = {
    'TIP3P': {'sigma_O': 3.1507, 'epsilon_O': 0.1521, 'q_O': -0.834, 'q_H': 0.417},
    'TIP4P': {'sigma_O': 3.1540, 'epsilon_O': 0.1550, 'q_O': 0.000, 'q_H': 0.520},  # Charge on M site
    'SPC/E': {'sigma_O': 3.1656, 'epsilon_O': 0.1554, 'q_O': -0.8476, 'q_H': 0.4238},
    'OPC':   {'sigma_O': 3.1668, 'epsilon_O': 0.2128, 'q_O': -0.8952, 'q_H': 0.4476},
}

r_water = np.linspace(2.0, 6, 300)

plt.figure(figsize=(10, 6))

for model, params in water_models.items():
    if model == 'TIP4P':
        continue  # Skip TIP4P for clarity (virtual site complicates comparison)
    
    sigma = params['sigma_O']
    epsilon = params['epsilon_O']
    U_LJ = 4 * epsilon * ((sigma/r_water)**12 - (sigma/r_water)**6)
    plt.plot(r_water, U_LJ, linewidth=2.5, label=model)

plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
plt.axvline(x=2.8, color='gray', linestyle=':', alpha=0.5, label='Typical O-O distance')
plt.xlabel('O-O Distance (Å)', fontsize=12)
plt.ylabel('LJ Energy (kcal/mol)', fontsize=12)
plt.title('Water Model Comparison: O-O LJ Interaction', fontsize=13, fontweight='bold')
plt.xlim(2.0, 6)
plt.ylim(-0.25, 0.5)
plt.legend(fontsize=11)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('water_models.png', dpi=150, bbox_inches='tight')
plt.show()

print("\n" + "="*60)
print("WATER MODEL PARAMETERS:")
print("="*60)
print(f"{'Model':<8} {'σ_O (Å)':<10} {'ε_O (kcal/mol)':<15} {'q_O (e)':<10} {'q_H (e)'}")
print("-" * 60)
for model, params in water_models.items():
    print(f"{model:<8} {params['sigma_O']:<10.4f} {params['epsilon_O']:<15.4f} "
          f"{params['q_O']:<10.4f} {params['q_H']:.4f}")
print("="*60)
print("\nExperimental targets at 25°C:")
print("  • Density: 0.997 g/cm³")
print("  • Heat of vaporization: 10.5 kcal/mol")
print("  • Dielectric constant: 78.4")
print("\n→ OPC is current best all-around performer")
print("→ TIP3P still widely used (fast, compatible)")
print("="*60)
```

### 2.6 When to Use Which Force Field

| System Type | Recommended Force Field |
|-------------|-------------------------|
| Proteins in water | AMBER ff19SB or CHARMM36m |
| RNA/DNA | AMBER OL15 or CHARMM36 |
| Lipid bilayers | CHARMM36 or Slipids |
| Carbohydrates | CHARMM36 or GLYCAM06 |
| Small organic molecules | OPLS-AA or GAFF (AMBER) |
| Polymers | OPLS-AA or TraPPE |
| Coarse-grained proteins | MARTINI |

**Golden Rule:** Use the force field that was validated for your specific system. Don't mix force fields (e.g., AMBER protein with CHARMM lipids) without careful validation.

---

## 3. Many-Body Potentials: The Embedded Atom Method (EAM)

### 3.1 Why Pair Potentials Fail for Metals

Consider copper metal. If we try to use LJ:
- LJ has a single minimum at $r_{min}$
- But copper atoms in FCC crystal have 12 nearest neighbors at one distance
- AND 6 second-nearest neighbors at a different distance
- A pair potential cannot simultaneously satisfy both

**The Fundamental Problem:**
In metals, the bonding energy of an atom depends on its **local environment** (how many neighbors it has), not just pairwise distances. This is a **many-body effect**.

**Experimental Evidence:**
- Vacancy formation energy (removing one atom) is NOT equal to breaking 12 bonds
- Surface atoms bind differently than bulk atoms
- Elastic constants cannot be fit with pair potentials

### 3.2 The EAM Concept: Embedding Energy

The Embedded Atom Method, developed by Murray Daw and Mike Baskes (1984), treats each atom as embedded in the electron density created by all other atoms.

**Total Energy:**
$$E_{total} = \sum_i F_i(\rho_i) + \frac{1}{2}\sum_{i \neq j} \phi_{ij}(r_{ij})$$

**Three functions define the potential:**

1. **Pair interaction** $\phi_{ij}(r)$: Repulsion between nuclei (like LJ repulsion)

2. **Electron density** $\rho_i$: Total density at atom i from all neighbors
   $$\rho_i = \sum_{j \neq i} f_j(r_{ij})$$
   where $f_j(r)$ is the electron density contribution from atom j

3. **Embedding energy** $F(\rho)$: Energy cost/gain of embedding atom i in density $\rho_i$

**Physical Interpretation:**
- High electron density → stronger bonding → lower energy
- Atoms with many neighbors (bulk) have higher $\rho$ than surface atoms
- $F(\rho)$ is fitted to reproduce cohesive energy, lattice constant, elastic constants

### 3.3 Interactive Example: EAM vs. Pair Potentials

```python
# Simplified EAM-like behavior visualization
r = np.linspace(2.0, 6.0, 300)

# Mock pair potential (repulsive part only)
phi_pair = 100 * np.exp(-2*(r - 2.5))

# Mock electron density contribution (decays with distance)
f_rho = np.exp(-(r - 2.5)**2 / 0.5)

# Mock embedding function (attractive, saturates at high density)
# For visualization: F(ρ) ≈ -A√ρ
rho_values = np.linspace(0, 5, 100)
F_embed = -2.5 * np.sqrt(rho_values)

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 5))

# --- Plot 1: Pair interaction ---
ax1.plot(r, phi_pair, 'b-', linewidth=2.5)
ax1.set_xlabel('Distance r (Å)', fontsize=12)
ax1.set_ylabel('φ(r) (eV)', fontsize=12)
ax1.set_title('Pair Repulsion φ(r)', fontsize=13, fontweight='bold')
ax1.grid(alpha=0.3)
ax1.set_ylim(-5, 50)

# --- Plot 2: Electron density ---
ax2.plot(r, f_rho, 'g-', linewidth=2.5)
ax2.set_xlabel('Distance r (Å)', fontsize=12)
ax2.set_ylabel('f(r)', fontsize=12)
ax2.set_title('Electron Density Contribution f(r)', fontsize=13, fontweight='bold')
ax2.grid(alpha=0.3)
ax2.set_ylim(0, 1.2)

# --- Plot 3: Embedding function ---
ax3.plot(rho_values, F_embed, 'r-', linewidth=2.5)
ax3.set_xlabel('Local Density ρ', fontsize=12)
ax3.set_ylabel('F(ρ) (eV)', fontsize=12)
ax3.set_title('Embedding Energy F(ρ)', fontsize=13, fontweight='bold')
ax3.grid(alpha=0.3)
ax3.axhline(y=0, color='k', linestyle='--', alpha=0.3)

plt.tight_layout()
plt.savefig('eam_components.png', dpi=150, bbox_inches='tight')
plt.show()

print("\n" + "="*60)
print("EAM PHYSICAL INTERPRETATION:")
print("="*60)
print("Bulk atom (12 neighbors):")
print("  ρ_bulk = 12 × f(r_nn) = HIGH")
print("  F(ρ_bulk) = VERY NEGATIVE (strong bonding)")
print("\nSurface atom (9 neighbors):")
print("  ρ_surface = 9 × f(r_nn) = LOWER")
print("  F(ρ_surface) = LESS NEGATIVE (weaker bonding)")
print("\n→ Surface atoms naturally have higher energy")
print("→ Explains surface tension, vacancy formation")
print("→ Cannot be captured by pair potentials!")
print("="*60)
```

### 3.4 EAM Parameter Files

EAM potentials are distributed as tabulated files (not analytical functions):

**File contains:**
- $\phi(r)$: pair interaction (tabulated values)
- $f(r)$: electron density function (tabulated)
- $F(\rho)$: embedding energy (tabulated)

**Standard formats:**
- **LAMMPS**: `.eam.alloy`, `.eam.fs` formats
- **DYNAMO**: Original format from Sandia National Labs

**Available elements:**
- FCC metals: Cu, Ag, Au, Ni, Pd, Pt, Al, Pb
- BCC metals: Fe, Cr, Mo, W, Ta, V
- HCP metals: Mg, Ti, Zr, Co
- Alloys: Many binary and ternary systems

```python
# Visualization: Why surface energy emerges from EAM
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(10, 6))

# Schematic: Atom coordination vs. energy
# Bulk (coordination 12), Surface (coordination 9), Adatom (coordination 1)
coordination = np.array([1, 3, 6, 9, 12])
# Mock electron density (proportional to coordination)
rho = coordination * 0.5
# Mock embedding energy F(ρ) ~ -A√ρ
F_energy = -3.0 * np.sqrt(rho)

ax.plot(coordination, F_energy, 'ro-', linewidth=2.5, markersize=12)
ax.axhline(y=F_energy[-1], color='blue', linestyle='--', alpha=0.5, 
           label=f'Bulk energy ({F_energy[-1]:.2f} eV)')
ax.axhline(y=F_energy[3], color='green', linestyle='--', alpha=0.5,
           label=f'Surface energy ({F_energy[3]:.2f} eV)')

# Annotate points
ax.annotate('Adatom', xy=(1, F_energy[0]), xytext=(2, F_energy[0]+0.5),
            fontsize=11, arrowprops=dict(arrowstyle='->', color='red'))
ax.annotate('Surface', xy=(9, F_energy[3]), xytext=(7, F_energy[3]-0.8),
            fontsize=11, arrowprops=dict(arrowstyle='->', color='green'))
ax.annotate('Bulk (FCC)', xy=(12, F_energy[4]), xytext=(10, F_energy[4]-0.8),
            fontsize=11, arrowprops=dict(arrowstyle='->', color='blue'))

ax.set_xlabel('Coordination Number', fontsize=12)
ax.set_ylabel('Embedding Energy F(ρ) (eV)', fontsize=12)
ax.set_title('EAM: Energy vs. Local Environment', fontsize=14, fontweight='bold')
ax.set_xlim(0, 13)
ax.set_ylim(-11, -2)
ax.legend(fontsize=11)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('eam_coordination.png', dpi=150, bbox_inches='tight')
plt.show()

print("\n" + "="*60)
print("SURFACE ENERGY FROM EAM:")
print("="*60)
print(f"Bulk atom (coord. 12): {F_energy[4]:.2f} eV")
print(f"Surface atom (coord. 9): {F_energy[3]:.2f} eV")
print(f"Energy difference: {F_energy[3] - F_energy[4]:.2f} eV/atom")
print("\n→ Surface atoms cost extra energy")
print("→ System minimizes surface area (surface tension)")
print("→ Realistic surface reconstruction emerges naturally")
print("="*60)
```

### 3.5 Applications of EAM

**What EAM Does Well:**
- ✅ Bulk crystal properties (lattice constant, cohesive energy)
- ✅ Elastic constants (C₁₁, C₁₂, C₄₄)
- ✅ Defects (vacancies, interstitials, stacking faults)
- ✅ Surfaces and interfaces
- ✅ Plastic deformation (dislocations, grain boundaries)
- ✅ Nanoparticles (size-dependent melting)

**What EAM Cannot Do:**
- ❌ Directional bonding (semiconductors like Si)
- ❌ Magnetic properties (requires spin)
- ❌ Chemical reactions (bonds don't break/form)
- ❌ Electronic structure (no band gaps)

**Typical Systems:**
- Metal nanoparticles (catalysis)
- Crack propagation in alloys
- Grain boundary sliding
- Radiation damage in reactor materials

---

## 4. Reactive Potentials: The ReaxFF Framework

### 4.1 The Challenge of Reactive Chemistry

Classical force fields have **fixed bonds**:
- C-C bond between atoms 1 and 2 is always present
- Cannot simulate bond breaking (combustion, pyrolysis)
- Cannot simulate bond formation (polymerization, catalysis)

**Reactive potentials** solve this by making bonds **dynamic**:
- Bond strength depends on instantaneous atomic positions
- Bonds can smoothly break and form during simulation
- No pre-defined reaction pathways needed

### 4.2 The Bond-Order Concept

**Key Idea:** Bond strength is not constant—it depends on the **bond order (BO)**, which changes continuously during the simulation.

**Definition:**
$$BO_{ij} = BO_{ij}^\sigma + BO_{ij}^\pi + BO_{ij}^{\pi\pi}$$

Each component depends on **interatomic distance**:

$$BO_{ij}^\sigma = \exp\left[p_{bo,1} \left(\frac{r_{ij}}{r_0^\sigma}\right)^{p_{bo,2}}\right]$$

**Physical Meaning:**
- $r_{ij}$ small (atoms close): $BO \approx 1$ (single bond), 2 (double), or 3 (triple)
- $r_{ij}$ large (atoms far): $BO \to 0$ (no bond)
- **Smooth transition** — bond smoothly weakens as atoms separate

```python
# Visualize bond order as a function of distance
r = np.linspace(1.0, 3.5, 300)

# Mock ReaxFF bond order parameters (C-C bond)
r0_sigma = 1.39  # Equilibrium bond length (Å)
p_bo1 = -0.1509
p_bo2 = 6.3890

# Bond order components
BO_sigma = np.exp(p_bo1 * (r / r0_sigma)**p_bo2)
BO_pi = np.exp(p_bo1 * (r / (r0_sigma*1.1))**p_bo2) * 0.7  # Mock π contribution
BO_total = BO_sigma + BO_pi

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# --- Plot 1: Bond order components ---
ax1.plot(r, BO_sigma, 'b-', linewidth=2.5, label='σ bond')
ax1.plot(r, BO_pi, 'g-', linewidth=2.5, label='π bond')
ax1.plot(r, BO_total, 'r-', linewidth=3, label='Total BO')

ax1.axhline(y=1, color='gray', linestyle='--', alpha=0.5, label='Single bond')
ax1.axhline(y=2, color='gray', linestyle=':', alpha=0.5, label='Double bond')
ax1.axvline(x=r0_sigma, color='orange', linestyle='--', alpha=0.5, label='r₀')

ax1.set_xlabel('C-C Distance (Å)', fontsize=12)
ax1.set_ylabel('Bond Order', fontsize=12)
ax1.set_title('ReaxFF: Bond Order vs. Distance', fontsize=13, fontweight='bold')
ax1.set_xlim(1.0, 3.5)
ax1.set_ylim(0, 2.5)
ax1.legend(fontsize=10)
ax1.grid(alpha=0.3)

# --- Plot 2: Bond energy from bond order ---
# Bond energy is a function of BO
De = 145  # kcal/mol (C-C bond dissociation energy)
E_bond = -De * BO_total * np.exp(-1.5 * (BO_total - 1))  # Mock functional form

ax2.plot(r, E_bond, 'purple', linewidth=3)
ax2.axhline(y=0, color='k', linestyle='--', alpha=0.3)
ax2.axvline(x=r0_sigma, color='orange', linestyle='--', alpha=0.5, label='r₀')

ax2.set_xlabel('C-C Distance (Å)', fontsize=12)
ax2.set_ylabel('Bond Energy (kcal/mol)', fontsize=12)
ax2.set_title('Bond Energy from Bond Order', fontsize=13, fontweight='bold')
ax2.set_xlim(1.0, 3.5)
ax2.set_ylim(-160, 10)
ax2.legend(fontsize=10)
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('reaxff_bond_order.png', dpi=150, bbox_inches='tight')
plt.show()

print("\n" + "="*60)
print("BOND ORDER INTERPRETATION:")
print("="*60)
print(f"At r = {r0_sigma} Å (equilibrium):")
print(f"  BO_total ≈ {BO_total[np.argmin(np.abs(r - r0_sigma))]:.2f}")
print(f"  → Strong C-C bond (single)")
print(f"\nAt r = 2.5 Å (stretched):")
print(f"  BO_total ≈ {BO_total[np.argmin(np.abs(r - 2.5))]:.2f}")
print(f"  → Weak partial bond")
print(f"\nAt r = 3.5 Å (dissociated):")
print(f"  BO_total ≈ {BO_total[-1]:.2f}")
print(f"  → No bond (non-bonded interaction takes over)")
print("\n→ Bond smoothly transitions from bonded to non-bonded")
print("→ No discontinuity — energies and forces are continuous")
print("="*60)
```

### 4.3 ReaxFF Energy Terms

The total ReaxFF energy contains many terms (30+ in full implementation):

$$E_{ReaxFF} = E_{bond} + E_{over} + E_{under} + E_{lp} + E_{val} + E_{tor} + E_{vdWaals} + E_{Coulomb}$$

**Key Components:**

1. **Bond energy** $E_{bond}$: Function of bond order
   $$E_{bond} = -D_e \cdot BO \cdot \exp[p_{be,1}(1-BO^{p_{be,2}})]$$

2. **Over-coordination penalty** $E_{over}$: Penalizes atoms with too many bonds
   - Prevents unphysical hyper-valent states
   
3. **Under-coordination stabilization** $E_{under}$: Stabilizes lone pairs, radicals
   
4. **Lone pair energy** $E_{lp}$: Accounts for electron lone pairs (e.g., on O, N)

5. **Valence angle** $E_{val}$: Depends on bond orders (stiff when BO is high)

6. **Torsion** $E_{tor}$: Dihedral barriers (depends on BO of central bond)

7. **van der Waals** $E_{vdWaals}$: Modified LJ-like term for non-bonded

8. **Coulomb** $E_{Coulomb}$: Electrostatics with **dynamic charges**

### 4.4 Charge Equilibration (QEq)

**The Problem:** In classical force fields, charges are fixed. But in reality:
- C=O carbon is more positive than C-C carbon
- Charges should change as bonds break/form

**ReaxFF Solution: Electronegativity Equalization**

At every MD step, atomic charges are recalculated by minimizing the **electrostatic energy**:

$$E_{elec} = \sum_i \left(\chi_i^0 q_i + \frac{1}{2}\eta_i q_i^2\right) + \sum_{i<j} \frac{q_i q_j}{r_{ij}}$$

Subject to: $\sum_i q_i = Q_{total}$ (charge conservation)

**Parameters:**
- $\chi_i^0$ = electronegativity of atom i
- $\eta_i$ = hardness (resistance to charge transfer)

**Result:** Charges redistribute to minimize energy, similar to how electrons would distribute in a real molecule.

```python
# Simplified QEq example: charge distribution
# Two atoms with different electronegativities

chi_C = 5.343   # Carbon electronegativity (eV)
chi_O = 8.741   # Oxygen electronegativity (eV)
eta_C = 5.063   # Carbon hardness (eV)
eta_O = 6.899   # Oxygen hardness (eV)

# Distance
r_CO = 1.43  # Angstroms

# Solve for charges (simplified - actual ReaxFF uses matrix inversion)
# Constraint: q_C + q_O = 0 (neutral molecule)
# Minimize: chi_C*q_C + 0.5*eta_C*q_C^2 + chi_O*q_O + 0.5*eta_O*q_O^2 + q_C*q_O/(r_CO)

# Analytical solution for 2-atom case:
k_e_eV = 14.4  # eV·Å/e² (Coulomb constant in eV units)
denominator = eta_C + eta_O + 2*k_e_eV/r_CO
q_C = (chi_O - chi_C) / denominator
q_O = -q_C

print("\n" + "="*60)
print("CHARGE EQUILIBRATION EXAMPLE: C=O BOND")
print("="*60)
print(f"Carbon electronegativity (χ_C): {chi_C:.3f} eV")
print(f"Oxygen electronegativity (χ_O): {chi_O:.3f} eV")
print(f"Distance (r_CO): {r_CO} Å")
print(f"\nEquilibrated charges:")
print(f"  q_C = {q_C:+.3f} e")
print(f"  q_O = {q_O:+.3f} e")
print(f"  Sum = {q_C + q_O:.6f} e (conserved)")
print("\n→ Oxygen is more electronegative → pulls charge")
print("→ C becomes positive, O becomes negative")
print("→ Charges adjust automatically as bond stretches/breaks")
print("="*60)

# Visualize how charges change with distance
r_range = np.linspace(1.2, 4.0, 100)
q_C_array = []
for r in r_range:
    denom = eta_C + eta_O + 2*k_e_eV/r
    q_C_array.append((chi_O - chi_C) / denom)

q_C_array = np.array(q_C_array)
q_O_array = -q_C_array

plt.figure(figsize=(10, 6))
plt.plot(r_range, q_C_array, 'b-', linewidth=2.5, label='q_C (carbon)')
plt.plot(r_range, q_O_array, 'r-', linewidth=2.5, label='q_O (oxygen)')
plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
plt.axvline(x=1.43, color='gray', linestyle=':', alpha=0.5, label='Equilibrium r_CO')

plt.xlabel('C-O Distance (Å)', fontsize=12)
plt.ylabel('Partial Charge (e)', fontsize=12)
plt.title('ReaxFF: Charge Redistribution During Bond Breaking', fontsize=13, fontweight='bold')
plt.legend(fontsize=11)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('qeq_example.png', dpi=150, bbox_inches='tight')
plt.show()

print("\nObservation:")
print(f"  At r = {r_range[0]:.1f} Å: q_C = {q_C_array[0]:+.3f}, q_O = {q_O_array[0]:+.3f}")
print(f"  At r = {r_range[-1]:.1f} Å: q_C = {q_C_array[-1]:+.3f}, q_O = {q_O_array[-1]:+.3f}")
print("  → As atoms separate, charges go to zero (neutral atoms)")
```

### 4.5 ReaxFF Parameter Files

ReaxFF force fields are **system-specific** and require extensive training:

**Training Data:**
- Quantum chemistry calculations (DFT energy landscapes)
- Bond dissociation energies
- Angle strain energies
- Reaction barriers
- Heat of formation

**File Format:** Plain text `.ff` file containing:
- ~50-100 parameters per element
- ~30 parameters per bond type
- Combining rules

**Example systems with ReaxFF:**
- Hydrocarbons (CHO): Combustion, pyrolysis
- Nitramines (CHNO): Explosives, propellants
- Silicon/Silicon oxide (SiO): Oxidation, glass
- Transition metals (NiFeCr): Catalysis, corrosion
- Polymers (polyethylene, PMMA): Thermal degradation

**CRITICAL:** ReaxFF is only as good as its training set. Using a hydrocarbon force field on a metal system will give garbage results.

### 4.6 Computational Cost

ReaxFF is **significantly more expensive** than classical force fields:

| Method | Scaling | Relative Cost | Timestep |
|--------|---------|---------------|----------|
| Classical FF | O(N) | 1× | 1-2 fs |
| ReaxFF | O(N) - O(N²) | 10-100× | 0.1-0.25 fs |
| DFT (VASP) | O(N³) | 10,000× | 0.5-1 fs |

**Why so slow?**
- Bond order calculation for ALL atom pairs (not just bonded)
- QEq charge calculation (matrix inversion) every step
- Shorter timestep needed (stiff bond-breaking forces)

**Practical Limits:**
- ~10,000 atoms for reasonable simulation times
- Nanoseconds timescale achievable
- Much faster than ab initio MD (femtoseconds only)

### 4.7 Applications and Success Stories

**ReaxFF Excels At:**
- ✅ Combustion mechanisms (hydrocarbon oxidation)
- ✅ Thermal decomposition of polymers
- ✅ Catalytic reactions on metal surfaces
- ✅ Battery electrolyte decomposition (SEI formation)
- ✅ Explosives (shock-induced chemistry)
- ✅ Reactive sputtering, thin film growth

**Example: Methane Combustion**

ReaxFF can simulate the full reaction network:
```
CH₄ + O₂ → CH₃ + OH
CH₃ + O₂ → CH₂O + OH
CH₂O → CO + H₂
CO + OH → CO₂ + H
H + O₂ → OH + O
...
```

All these reactions emerge automatically—no pre-defined pathways. The simulation "discovers" the reaction mechanism.

```python
# Visualization: Schematic reaction pathway in ReaxFF
import matplotlib.pyplot as plt
import numpy as np

# Mock reaction coordinate
reaction_coord = np.linspace(0, 10, 300)

# Mock potential energy surface with transition states
# Reactants → TS1 → Intermediate → TS2 → Products
E_surface = (
    50 * np.exp(-((reaction_coord - 0.5)/0.5)**2) +  # Reactant well
    100 * np.exp(-((reaction_coord - 3)/0.3)**2) +   # TS1
    30 * np.exp(-((reaction_coord - 5)/0.7)**2) +    # Intermediate
    90 * np.exp(-((reaction_coord - 7.5)/0.3)**2) -  # TS2
    80  # Baseline shift to products
)

plt.figure(figsize=(10, 6))
plt.plot(reaction_coord, E_surface, 'b-', linewidth=3)

# Annotate key points
plt.annotate('Reactants\n(CH₄ + O₂)', xy=(0.5, 50), xytext=(0.5, 70),
            ha='center', fontsize=11, fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
plt.annotate('TS1', xy=(3, 100), xytext=(3, 120),
            ha='center', fontsize=11, fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='red', lw=1.5))
plt.annotate('Intermediate\n(CH₃ + OH)', xy=(5, 30), xytext=(5, -10),
            ha='center', fontsize=11, fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
plt.annotate('Products\n(CO₂ + H₂O)', xy=(9.5, -30), xytext=(9.5, -50),
            ha='center', fontsize=11, fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='green', lw=1.5))

plt.xlabel('Reaction Coordinate', fontsize=12)
plt.ylabel('Potential Energy (kcal/mol)', fontsize=12)
plt.title('ReaxFF: Reactive Potential Energy Surface', fontsize=14, fontweight='bold')
plt.ylim(-60, 130)
plt.grid(alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig('reaxff_reaction_pathway.png', dpi=150, bbox_inches='tight')
plt.show()

print("\n" + "="*60)
print("REAXFF REACTIVE SIMULATION:")
print("="*60)
print("What happens during simulation:")
print("  1. Atoms start as CH₄ and O₂ molecules (bonded)")
print("  2. Thermal fluctuations bring them together")
print("  3. C-H bond weakens (BO decreases)")
print("  4. O-H bond forms (BO increases)")
print("  5. System passes through transition state")
print("  6. New molecules (CO₂, H₂O) emerge")
print("\n→ No reaction pathway pre-programmed")
print("→ Chemistry emerges from bond order dynamics")
print("→ Can discover unexpected reaction mechanisms")
print("="*60)
```

---

## 5. Comparison and Summary

### 5.1 Choosing the Right Potential

| System | Classical FF | EAM | ReaxFF |
|--------|--------------|-----|--------|
| Proteins, DNA | ✅ AMBER/CHARMM | ❌ | ❌ |
| Small molecules in solvent | ✅ OPLS | ❌ | ❌ |
| Metal nanoparticles | ❌ | ✅ | ❌ |
| Metal alloys, defects | ❌ | ✅ | ❌ |
| Combustion | ❌ | ❌ | ✅ |
| Polymer pyrolysis | ❌ | ❌ | ✅ |
| Catalytic reactions | ❌ | ❌ | ✅ |
| Battery SEI formation | ❌ | ❌ | ✅ |

### 5.2 The Computational Cost Ladder

```
Accuracy vs. Speed Trade-off:

Classical FF (fastest)
    ↓  10× slower
EAM (metals only)
    ↓  10-100× slower  
ReaxFF (reactive)
    ↓  100× slower
Ab Initio MD (DFT)
    ↓  1000× slower
Quantum Chemistry (CCSD(T))
```

**The Golden Rule:** Use the simplest potential that captures your physics.

### 5.3 Validation is Everything

Regardless of which potential you choose:

**Always validate against experiments:**
- Structure (X-ray, neutron scattering)
- Thermodynamics (density, heat capacity, phase diagram)
- Dynamics (diffusion, viscosity, relaxation times)

**Never blindly trust simulation results.** A beautiful trajectory from a bad force field is worthless.

---

## 6. Practical Resources

### 6.1 Where to Get Parameters

**Classical Force Fields:**
- AMBER: [AmberTools](http://ambermd.org/AmberTools.php) (free)
- CHARMM: [CHARMM-GUI](https://www.charmm-gui.org/) (free, excellent web interface)
- OPLS: [GROMACS force field repository](http://www.gromacs.org/)

**EAM Potentials:**
- NIST repository: [https://www.ctcms.nist.gov/potentials/](https://www.ctcms.nist.gov/potentials/)
- LAMMPS distribution includes many EAM files
- OpenKIM: [https://openkim.org/](https://openkim.org/) (standardized database)

**ReaxFF:**
- SCM (Software for Chemistry & Materials): Commercial, curated force fields
- Literature: Many published parameter sets
- LAMMPS `examples/reax` directory

### 6.2 Software Support

| Potential Type | LAMMPS | GROMACS | AMBER | NAMD |
|----------------|--------|---------|-------|------|
| Classical FF | ✅ | ✅ | ✅ | ✅ |
| Water models | ✅ | ✅ | ✅ | ✅ |
| EAM | ✅ | ❌ | ❌ | ❌ |
| ReaxFF | ✅ | ❌ | ❌ | ❌ |

**LAMMPS is the most versatile** for materials science (supports all potential types).

---

## 7. Key Takeaways

**What You Should Remember:**

1. **Force fields are complete parameter sets** - Don't mix AMBER and CHARMM!

2. **Atom types are chemical environments** - CT ≠ CA ≠ C

3. **Water models matter** - TIP3P vs. TIP4P can change your results

4. **EAM captures many-body effects** - Essential for metals, cannot use pair potentials

5. **ReaxFF enables reactive simulations** - Bond order + QEq allows chemistry

6. **Computational cost varies 100×** - Classical << EAM < ReaxFF << DFT

7. **Validation is non-negotiable** - Beautiful simulations mean nothing without experimental validation

**In the tutorials, you'll:**
- Use LAMMPS to run simulations with these potentials
- See concrete examples of what each potential can and cannot predict
- Learn to read parameter files and troubleshoot errors
- Compare results across different potential types

---

## 8. Looking Ahead: Machine Learning Potentials

The future of MD potentials is **machine learning**:

**Modern Approaches:**
- Neural Network Potentials (NNP)
- Gaussian Approximation Potentials (GAP)
- Graph Neural Networks (SchNet, DimeNet++)

**The Promise:**
- Near-DFT accuracy
- Classical MD speed (1000× faster than DFT)
- Transferable across systems

**The Challenge:**
- Require massive training datasets
- Extrapolation beyond training data is unreliable
- Still under active development

This is an active research area—stay tuned for the next decade of advances!