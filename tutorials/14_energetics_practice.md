# Tutorial 14: Energetics in Practice

**Objective:** Compute energetic quantities using ASE and ML potentials. Move from 0 K energy differences to finite-temperature free energies.

:::{note} Prerequisites
- Tutorial 13 (MACE basics and ASE interface)
- Lecture 11 (energetics and free energy concepts)
- The `molsimclass` conda environment on Anvil
:::

---

## Setup

### Working Directory

Open a terminal in your Jupyter session and create a working directory:

```bash
cd $SCRATCH
mkdir -p tutorial14
cd tutorial14
```

### Calculators Used in This Tutorial

This tutorial uses two ML calculators:

- **MACE-MP-0** (`mace_mp`) — for bulk metals. Already familiar from Tutorial 13.
- **MACE-OFF** (`mace_off`) — for organic molecules. Same package, different model weights. Downloads automatically on first use (~150 MB).

Run this cell first to verify both are available:

```python
import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.build import bulk, molecule
from mace.calculators import mace_mp, mace_off

print("All imports successful.")
```

---

## Part 1: Geometry Optimization and Cohesive Energy

**System:** Bulk FCC copper using MACE-MP-0.

**Concepts:** Energy minimization in practice, cohesive energy as an energy difference.

### Cell 1: Load Calculator and Build Structure

```python
# Load MACE-MP-0 (medium accuracy, float64 for better precision)
calc = mace_mp(model="medium", dispersion=False, default_dtype="float64")
print("MACE-MP-0 loaded.")

# Build FCC copper at the experimental lattice constant
cu = bulk('Cu', crystalstructure='fcc', a=3.615, cubic=True)
cu.calc = calc
print(f"Cu unit cell: {len(cu)} atoms")
```

### Cell 2: Forces Before and After Optimization

Before minimization, a perfect crystal should already have near-zero forces. After a small perturbation, they will not be.

```python
# Perfect crystal: forces should be ~0
energy_before = cu.get_potential_energy()
forces_before  = cu.get_forces()
print(f"Energy (perfect crystal):        {energy_before:.6f} eV")
print(f"Max force (perfect crystal):     {np.max(np.abs(forces_before)):.2e} eV/Å")
print()

# Perturb positions slightly to create non-zero forces
cu_perturbed = cu.copy()
cu_perturbed.rattle(stdev=0.05, seed=42)
cu_perturbed.calc = calc

forces_perturbed = cu_perturbed.get_forces()
print(f"Max force (after rattle 0.05 Å): {np.max(np.abs(forces_perturbed)):.4f} eV/Å")
print("→ Optimization is needed to restore the minimum.")
```

### Cell 3: Run Geometry Optimization

```python
from ase.optimize import BFGS

cu_opt = cu_perturbed.copy()
cu_opt.calc = calc

opt = BFGS(cu_opt, logfile=None)
opt.run(fmax=0.005)

print(f"Converged in {opt.nsteps} steps.")
print(f"Max force after optimization: {np.max(np.abs(cu_opt.get_forces())):.2e} eV/Å")
print(f"Energy: {cu_opt.get_potential_energy():.6f} eV")
```

### Cell 4: Cohesive Energy

Cohesive energy requires the energy of an isolated Cu atom. We place a single atom in a large box so no periodic images fall within the potential cutoff.

```python
# Isolated Cu atom in a large vacuum box
cu_atom = Atoms('Cu', positions=[(0, 0, 0)], cell=[20, 20, 20], pbc=True)
cu_atom.calc = calc
E_atom = cu_atom.get_potential_energy()

n_atoms    = len(cu_opt)
E_crystal  = cu_opt.get_potential_energy()
E_per_atom = E_crystal / n_atoms

# Cohesive energy: energy released per atom on condensation
E_coh = E_atom - E_per_atom    # positive = crystal is more stable than isolated atoms

print(f"E(isolated Cu atom):  {E_atom:.6f} eV")
print(f"E(crystal)/atom:      {E_per_atom:.6f} eV")
print()
print(f"Cohesive energy (MACE): {E_coh:.4f} eV/atom")
print(f"Experimental:            3.49  eV/atom")
print(f"Error: {(E_coh - 3.49)/3.49 * 100:+.1f}%")
```

---

## Part 2: Equation of State and Bulk Modulus

**Concept:** Scan energy as a function of volume, then fit the Birch-Murnaghan equation of state to extract the equilibrium lattice constant $a_0$, static energy $E_0$, and bulk modulus $B_0$.

The bulk modulus measures resistance to uniform compression:

$$B_0 = -V \frac{\partial P}{\partial V} \bigg|_{T=0}$$

It is related to the curvature of the energy–volume curve at the minimum — a narrow, steep well means a large $B_0$ and a stiff material.

### Cell 5: Volume Scan

```python
print("Scanning volumes... (takes ~2 minutes)")

a_values = np.linspace(3.40, 3.80, 20)
volumes  = []
energies = []

for a in a_values:
    cu_test = bulk('Cu', crystalstructure='fcc', a=a, cubic=True)
    cu_test.calc = calc
    volumes.append(cu_test.get_volume())
    energies.append(cu_test.get_potential_energy())
    print(f"  a = {a:.3f} Å  →  E = {energies[-1]:.5f} eV")

volumes  = np.array(volumes)
energies = np.array(energies)
n_cell   = len(bulk('Cu', crystalstructure='fcc', a=3.615, cubic=True))
```

### Cell 6: Birch-Murnaghan Fit

```python
from ase.eos import EquationOfState

eos = EquationOfState(volumes, energies, eos='birchmurnaghan')
v0, e0, B = eos.fit()

# Convert B from eV/Å³ to GPa (1 eV/Å³ = 160.2176 GPa)
B_GPa  = B * 160.2176
a0_fit = (v0 / n_cell * 4)**(1/3)   # back-calculate lattice constant for FCC

print("="*50)
print("BIRCH-MURNAGHAN FIT: Copper")
print("="*50)
print(f"Lattice constant a₀:  {a0_fit:.4f} Å    (exp: 3.615 Å)")
print(f"Cohesive energy E₀:   {e0/n_cell:.4f} eV/atom")
print(f"Bulk modulus B₀:      {B_GPa:.1f} GPa   (exp: 137 GPa)")
```

### Cell 7: EOS Plot

```python
from ase.eos import birchmurnaghan

fig, ax = plt.subplots(figsize=(9, 6))

# Raw data (per atom)
ax.scatter(volumes / n_cell, energies / n_cell,
           s=60, color='steelblue', zorder=5, label='MACE-MP-0')

# Fitted curve
v_fit = np.linspace(volumes.min(), volumes.max(), 200)
e_fit = birchmurnaghan(v_fit, e0, B, 4.0, v0)
ax.plot(v_fit / n_cell, e_fit / n_cell, 'b-', lw=2, label='Birch-Murnaghan fit')

ax.axvline(v0 / n_cell, color='red', ls='--', lw=1.5,
           label=f'V₀ = {v0/n_cell:.2f} Å³/atom')
ax.axhline(e0 / n_cell, color='gray', ls=':', lw=1)

ax.set_xlabel('Volume per atom (Å³)', fontsize=12)
ax.set_ylabel('Energy per atom (eV)', fontsize=12)
ax.set_title('Copper Equation of State', fontsize=13, fontweight='bold')
ax.legend(fontsize=11)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('cu_eos.png', dpi=150)
plt.show()

print(f"Summary: a₀ = {a0_fit:.4f} Å,  B₀ = {B_GPa:.1f} GPa")
```

---

## Part 3: Molecular Free Energy — H₂ Dissociation ΔG(T)

**Concept:** Use MACE-OFF for a molecular system. Compute the potential energy curve for H₂, then calculate $G(T)$ using ideal gas statistical mechanics to show that H₂ is thermodynamically stable at 300 K but dissociates spontaneously at extreme temperatures.

### Cell 8: Load MACE-OFF

MACE-OFF is trained on molecular data and is more appropriate than MACE-MP-0 for gas-phase molecules.

```python
print("Loading MACE-OFF (small model)...")
calc_mol = mace_off(model="small", default_dtype="float64")
print("MACE-OFF loaded.")
```

### Cell 9: Optimize H₂ Geometry

```python
from ase.optimize import BFGS

# Build H₂ at approximate experimental bond length
h2 = Atoms('H2', positions=[(0, 0, 0), (0, 0, 0.74)],
           cell=[15, 15, 15], pbc=False)
h2.calc = calc_mol

opt = BFGS(h2, logfile=None)
opt.run(fmax=0.001)

bond_length = h2.get_distance(0, 1)
E_h2        = h2.get_potential_energy()

print(f"Optimized H₂ bond length: {bond_length:.4f} Å   (exp: 0.741 Å)")
print(f"H₂ total energy:          {E_h2:.6f} eV")
```

### Cell 10: Potential Energy Curve

Scan the H–H distance to visualize the bond and estimate the dissociation energy.

```python
print("Scanning H–H potential energy curve...")

d_values = np.linspace(0.5, 4.0, 50)
E_curve  = []

for d in d_values:
    h2_test = Atoms('H2', positions=[(0, 0, 0), (0, 0, d)],
                    cell=[15, 15, 15], pbc=False)
    h2_test.calc = calc_mol
    E_curve.append(h2_test.get_potential_energy())

E_curve  = np.array(E_curve)
E_dissoc = E_curve[-1]          # large separation ≈ 2 × E(H atom)
D_e      = E_dissoc - E_curve.min()

fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(d_values, E_curve - E_dissoc, 'b-', lw=2.5, label='MACE-OFF')
ax.axhline(0, color='gray', ls='--', lw=1, label='Dissociation limit')
ax.axvline(bond_length, color='red', ls=':', lw=1.5,
           label=f'r₀ = {bond_length:.3f} Å')

ax.set_xlabel('H–H Distance (Å)', fontsize=12)
ax.set_ylabel('Relative Energy (eV)', fontsize=12)
ax.set_title('H₂ Potential Energy Curve (MACE-OFF)', fontsize=13, fontweight='bold')
ax.legend(fontsize=11)
ax.grid(alpha=0.3)
ax.set_xlim(0.5, 4.0)
ax.set_ylim(-5.5, 1.0)
plt.tight_layout()
plt.savefig('h2_pes.png', dpi=150)
plt.show()

print(f"\nClassical dissociation energy Dₑ: {D_e:.3f} eV/bond")
print(f"Experimental Dₑ:                  4.52 eV/bond")
```

### Cell 11: Vibrational Frequencies

`IdealGasThermo` needs the vibrational frequencies of H₂. ASE computes them by finite-difference differentiation of forces.

```python
from ase.vibrations import Vibrations

h2.calc = calc_mol   # re-attach calculator before running Vibrations

vib = Vibrations(h2, name='h2_vib', delta=0.01)
vib.run()

vib_energies_all = vib.get_energies()
print("All vibrational modes (eV):")
for i, e in enumerate(vib_energies_all):
    print(f"  Mode {i}: {e:.4f} eV  ({e.real * 8065.5:.0f} cm⁻¹)")

# Keep only real, non-zero modes (vibrational stretch; drop translations/rotations)
vib_energies = [e.real for e in vib_energies_all if e.real > 0.05]
print(f"\nRetained vibrational modes: {len(vib_energies)}")
print(f"H₂ stretch: {vib_energies[0]:.4f} eV  ({vib_energies[0]*8065.5:.0f} cm⁻¹)")
print(f"Experimental:               0.516 eV  (4161 cm⁻¹)")
```

### Cell 12: Gibbs Free Energy — IdealGasThermo

```python
from ase.thermochemistry import IdealGasThermo

P = 101325.0   # Standard pressure (Pa)

# --- H₂ molecule ---
# potentialenergy must be passed explicitly; IdealGasThermo defaults to 0 otherwise
thermo_h2 = IdealGasThermo(
    vib_energies    = vib_energies,
    geometry        = 'linear',
    atoms           = h2,
    symmetrynumber  = 2,          # homonuclear diatomic: σ = 2
    spin            = 0,          # singlet ground state
    potentialenergy = E_h2,
)

# --- Isolated H atom ---
h_atom = Atoms('H', positions=[(0, 0, 0)], cell=[15, 15, 15], pbc=False)
h_atom.calc = calc_mol
E_h_atom = h_atom.get_potential_energy()

thermo_h = IdealGasThermo(
    vib_energies    = [],           # monatomic: no vibrations
    geometry        = 'monatomic',
    atoms           = h_atom,
    symmetrynumber  = 1,
    spin            = 0.5,          # doublet: one unpaired electron
    potentialenergy = E_h_atom,
)

# --- ΔG(T) = 2G(H, T) − G(H₂, T) ---
T_values = np.linspace(300, 5000, 100)
delta_G  = []

for T in T_values:
    G_h2 = thermo_h2.get_gibbs_energy(T, P, verbose=False)
    G_h  = thermo_h.get_gibbs_energy( T, P, verbose=False)
    delta_G.append(2.0 * G_h - G_h2)

delta_G = np.array(delta_G)

# Crossover temperature
cross_idx = np.where(np.diff(np.sign(delta_G)))[0]
if len(cross_idx) > 0:
    i = cross_idx[0]
    T_star = T_values[i] + delta_G[i] * (T_values[i] - T_values[i+1]) \
             / (delta_G[i+1] - delta_G[i])
    print(f"ΔG = 0 at T* ≈ {T_star:.0f} K ({T_star - 273:.0f} °C)")
    print("Above this temperature, dissociation is thermodynamically spontaneous.")
else:
    T_star = None
    print("No crossover found below 5000 K.")
```

### Cell 13: Plot ΔG(T) and Interpret

```python
fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(T_values, delta_G, 'b-', lw=2.5,
        label=r'$\Delta G = 2G(\mathrm{H}) - G(\mathrm{H_2})$')
ax.axhline(0, color='k', ls='--', lw=1.5)

if T_star:
    ax.axvline(T_star, color='red', ls=':', lw=2,
               label=f'$T^*$ ≈ {T_star:.0f} K')

ax.fill_between(T_values, delta_G, 0,
                where=delta_G > 0, alpha=0.12, color='blue',
                label=r'$\Delta G > 0$: H₂ stable')
ax.fill_between(T_values, delta_G, 0,
                where=delta_G < 0, alpha=0.12, color='red',
                label=r'$\Delta G < 0$: dissociation spontaneous')

ax.annotate('300 K: H₂ very stable',
            xy=(300, delta_G[0]), xytext=(700, delta_G[0] * 0.75),
            fontsize=10, arrowprops=dict(arrowstyle='->', color='black'))

ax.set_xlabel('Temperature (K)', fontsize=12)
ax.set_ylabel('ΔG (eV per bond)', fontsize=12)
ax.set_title(r'H$_2$ Dissociation Free Energy vs Temperature',
             fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('h2_dissociation_deltaG.png', dpi=150)
plt.show()

print("""
INTERPRETATION
==============
At 300 K:  ΔG >> 0  →  H₂ is strongly stable.
           Breaking the bond costs ~4 eV; entropy gain is insufficient.

At T*:     ΔG = 0   →  Dissociation becomes spontaneous.
           (~3000–4000 K, plasma conditions)

Why does ΔS drive dissociation at high T?
  H₂ → 2H produces two independent particles instead of one.
  Each H atom gains full translational freedom.
  ΔS > 0, so −TΔS becomes large and negative at high T,
  eventually overcoming the positive ΔH.

General rule:  T* ≈ ΔH / ΔS
""")
```

---

## Part 4: LAMMPS Reference

The calculations above used ASE with MACE. The LAMMPS equivalents are shown here for reference — useful for running energy minimizations on large systems or as part of a workflow that feeds into MD.

::::{dropdown} LAMMPS: Geometry Optimization (minimize command)

**Input file: `minimize_cu.in`**

```lammps
# Copper cohesive energy via energy minimization
units       metal
atom_style  atomic
boundary    p p p

# Build FCC lattice
lattice     fcc 3.615
region      box block 0 3 0 3 0 3
create_box  1 box
create_atoms 1 box

# EAM potential
pair_style  eam
pair_coeff  1 1 /anvil/projects/x-chm250117/potentials/Cu_u3.eam

# Minimize with conjugate gradient
min_style   cg
minimize    1e-8 1e-10 10000 100000

# Report energy
variable E_total equal pe
variable N_atoms equal natoms
variable E_atom  equal v_E_total/v_N_atoms

print "Total energy:  ${E_total} eV"
print "Energy/atom:   ${E_atom} eV"
```

**Key LAMMPS minimize keywords:**

`min_style cg` — conjugate gradient (default, good for crystals)

`min_style fire` — FIRE algorithm (better for large disordered systems)

`minimize etol ftol maxiter maxeval`
- `etol`: relative energy change tolerance (dimensionless)
- `ftol`: force convergence (eV/Å in metal units)
- `maxiter`: maximum iterations
- `maxeval`: maximum force evaluations

**Pre-run outputs** for EAM Cu are in the class directory:
```
/anvil/projects/x-chm250117/class_examples/tutorial14/lammps_minimize/
```
::::

::::{dropdown} LAMMPS: Reading Minimization Output in Python

```python
import numpy as np
import matplotlib.pyplot as plt

# Parse LAMMPS log for potential energy during minimization
energies_lammps = []
with open('minimize_cu.out', 'r') as f:
    reading = False
    for line in f:
        if 'Step' in line and 'PotEng' in line:
            reading = True
            continue
        if reading:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    energies_lammps.append(float(parts[1]))
                except ValueError:
                    reading = False

energies_lammps = np.array(energies_lammps)

plt.figure(figsize=(8, 4))
plt.plot(energies_lammps, 'b-', lw=2)
plt.xlabel('Minimization Step')
plt.ylabel('Potential Energy (eV)')
plt.title('LAMMPS Minimization Convergence')
plt.grid(alpha=0.3)
plt.show()
```
::::

---

## Summary

| Part | System | Calculator | Key Result |
|:-----|:-------|:-----------|:-----------|
| 1 | Cu bulk | MACE-MP-0 | Cohesive energy vs experiment |
| 2 | Cu bulk | MACE-MP-0 | EOS: lattice constant, bulk modulus |
| 3 | H₂ molecule | MACE-OFF | ΔG(T): entropy drives dissociation at high T |

The central lesson: energy alone (Parts 1–2) answers "which structure is stable at 0 K?" Free energy (Part 3) answers "which state is stable at temperature T?" Both are necessary for a complete picture.

---

## Discussion Questions

1. In Part 1, you placed an isolated Cu atom in a 20 Å box to compute $E_{\text{atom}}$. Why does the box need to be large? What would happen if you used a 5 Å box?

2. In Part 2, the bulk modulus $B_0$ comes from the curvature of the energy–volume curve. Qualitatively, which would have a larger $B_0$: diamond or lead? Why?

3. In Part 3, the entropy change $\Delta S$ for H₂ dissociation is positive. Explain physically why two H atoms have more entropy than one H₂ molecule.

4. The static calculation gives $\Delta E \approx +4.5$ eV for H₂ dissociation but $T^*$ is only ~3000 K. Why does it take such a high temperature to reach $\Delta G = 0$, even though $\Delta S > 0$ and should help dissociation at any temperature above 0 K?
