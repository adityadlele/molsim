# Classical Interatomic Potentials: From Pair Interactions to Molecular Force Fields

---

## 1. Introduction: The Energy Landscape

### 1.1 The Central Problem of Molecular Dynamics

In the previous chapter, we established that Molecular Dynamics (MD) is fundamentally the integration of Newton's equations of motion:

$$\mathbf{F}_i = m_i \mathbf{a}_i = m_i \frac{d^2\mathbf{r}_i}{dt^2}$$

But to calculate the force, we need to know the **Potential Energy Surface** $U(\mathbf{r}^N)$, which describes how the energy of the system depends on the positions of all $N$ atoms. The force on atom $i$ is simply:

$$\mathbf{F}_i = -\nabla_i U(\mathbf{r}^N) = -\frac{\partial U}{\partial \mathbf{r}_i}$$

**The Challenge:**  
To compute $U(\mathbf{r}^N)$ from first principles (quantum mechanics), we would need to solve the Schrödinger equation for all electrons in the system at every timestep. For a protein with 10,000 atoms, this is computationally impossible with current technology.

**The Solution:**  
We use **MD potentials** (also called force fields) that approximate the potential energy using simple mathematical functions with parameters fitted to experimental data or high-level quantum calculations.

### 1.2 The Philosophy: A Hierarchy of Approximations

Different systems require different levels of sophistication. The art of molecular simulation is choosing the simplest potential that captures the physics you care about.

| System Type | Key Physics | Potential Complexity |
|-------------|-------------|---------------------|
| Noble gases (Ar, Ne) | Weak van der Waals forces | Pair potential (Lennard-Jones) |
| Simple molecules (N₂, O₂) | Covalent bonds | Pair + bond stretching |
| Water, small molecules | Bonds + angles | Pair + bonds + angles |
| Proteins, polymers | Full molecular structure | Bonds + angles + dihedrals + cross-terms |
| Reactive chemistry | Bond breaking/forming | Reactive potentials (ReaxFF) |
| Metals, alloys | Many-body interactions | Embedded Atom Method (EAM) |

**Our Approach:**  
We will build up complexity step by step, starting with the simplest case (two argon atoms) and progressing to full molecular force fields.

---

## 2. Pair Potentials: The Lennard-Jones Model

### 2.1 Physical Basis: Why Atoms Interact

Even chemically inert atoms like Argon experience forces when they approach each other. These forces arise from two quantum mechanical effects:

1. **Pauli Repulsion (Short Range):**  
   When electron clouds overlap, the Pauli Exclusion Principle creates a strong repulsive force. Energy rises steeply as atoms get too close.

2. **London Dispersion (Long Range):**  
   Instantaneous fluctuations in the electron cloud create temporary dipoles, which induce dipoles in neighboring atoms. This creates a weak attractive force proportional to $r^{-6}$.

### 2.2 The Lennard-Jones (12-6) Potential

The most famous empirical potential combines these two effects:

$$U_{LJ}(r) = 4\epsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6\right]$$

**Parameters:**
- $\epsilon$ = depth of the potential well (energy scale, typically in kcal/mol or kJ/mol)
- $\sigma$ = distance at which $U(r) = 0$ (length scale, typically in Angstroms)

**Key Distances:**
- $r_{min} = 2^{1/6}\sigma \approx 1.122\sigma$ : equilibrium distance where $U(r)$ is minimized
- $U(r_{min}) = -\epsilon$ : minimum energy (attraction)

**Why these exponents?**
- The $r^{-6}$ attraction is physically motivated (London dispersion).
- The $r^{-12}$ repulsion is chosen for computational convenience: $r^{12} = (r^6)^2$ can be calculated with one multiplication.

### 2.3 The Force from the Potential

The force between two atoms is the negative derivative:

$$F_{LJ}(r) = -\frac{dU}{dr} = \frac{24\epsilon}{\sigma}\left[2\left(\frac{\sigma}{r}\right)^{13} - \left(\frac{\sigma}{r}\right)^7\right]$$

**Key Features:**
- $F(r_{min}) = 0$ : no net force at equilibrium
- $F(r < r_{min}) > 0$ : repulsive (push apart)
- $F(r > r_{min}) < 0$ : attractive (pull together)
- $F(r \to \infty) \to 0$ : force vanishes at large separation

### 2.4 Interactive Example: Exploring the LJ Potential

Let's visualize how the parameters $\epsilon$ and $\sigma$ affect the potential energy curve.

```python
import numpy as np
import matplotlib.pyplot as plt

def lennard_jones(r, epsilon, sigma):
    """
    Compute Lennard-Jones potential energy.
    
    Parameters:
    -----------
    r : array
        Distance between atoms
    epsilon : float
        Depth of potential well (energy units)
    sigma : float
        Zero-crossing distance (length units)
    
    Returns:
    --------
    U : array
        Potential energy at each distance
    """
    return 4 * epsilon * ((sigma/r)**12 - (sigma/r)**6)

def lj_force(r, epsilon, sigma):
    """
    Compute Lennard-Jones force (magnitude).
    Positive = repulsive, Negative = attractive
    """
    return 24 * epsilon / sigma * (2*(sigma/r)**13 - (sigma/r)**7)

# Create distance array (avoid r=0)
r = np.linspace(0.8, 3.5, 500)

# Standard Argon parameters
epsilon_Ar = 1.0  # Normalized units
sigma_Ar = 1.0

# Calculate potential and force
U_standard = lennard_jones(r, epsilon_Ar, sigma_Ar)
F_standard = lj_force(r, epsilon_Ar, sigma_Ar)

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# --- Plot 1: Potential Energy ---
ax1.plot(r, U_standard, 'b-', linewidth=2.5, label='Standard LJ')
ax1.axhline(y=0, color='k', linestyle='--', alpha=0.3)
ax1.axhline(y=-epsilon_Ar, color='r', linestyle='--', alpha=0.5, 
            label=f'Well depth = -ε = {-epsilon_Ar}')
ax1.axvline(x=2**(1/6)*sigma_Ar, color='g', linestyle='--', alpha=0.5,
            label=f'r_min = 2^(1/6)σ ≈ {2**(1/6):.3f}')
ax1.axvline(x=sigma_Ar, color='orange', linestyle=':', alpha=0.5,
            label=f'σ = {sigma_Ar}')

ax1.set_xlim(0.8, 3.5)
ax1.set_ylim(-1.5, 2)
ax1.set_xlabel('Distance r/σ', fontsize=12)
ax1.set_ylabel('Energy U/ε', fontsize=12)
ax1.set_title('Lennard-Jones Potential Energy Curve', fontsize=14, fontweight='bold')
ax1.legend(fontsize=10)
ax1.grid(alpha=0.3)

# --- Plot 2: Force ---
ax2.plot(r, F_standard, 'r-', linewidth=2.5, label='F(r)')
ax2.axhline(y=0, color='k', linestyle='--', alpha=0.3)
ax2.axvline(x=2**(1/6)*sigma_Ar, color='g', linestyle='--', alpha=0.5,
            label='F=0 at r_min')
ax2.fill_between(r, 0, F_standard, where=(F_standard > 0), 
                  alpha=0.2, color='red', label='Repulsion')
ax2.fill_between(r, 0, F_standard, where=(F_standard < 0), 
                  alpha=0.2, color='blue', label='Attraction')

ax2.set_xlim(0.8, 3.5)
ax2.set_ylim(-2, 4)
ax2.set_xlabel('Distance r/σ', fontsize=12)
ax2.set_ylabel('Force F×σ/ε', fontsize=12)
ax2.set_title('Lennard-Jones Force Curve', fontsize=14, fontweight='bold')
ax2.legend(fontsize=10)
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('LJ_potential_and_force.png', dpi=150, bbox_inches='tight')
plt.show()

print("\n" + "="*60)
print("INTERPRETATION:")
print("="*60)
print(f"• At r < r_min ({2**(1/6):.3f}σ): Atoms repel (positive force)")
print(f"• At r = r_min ({2**(1/6):.3f}σ): Equilibrium (zero force)")
print(f"• At r > r_min: Atoms attract (negative force)")
print(f"• Well depth: {epsilon_Ar} ε (binding energy)")
print("="*60)
```

### 2.5 Effect of Parameters: Interactive Exploration

Now let's see how changing $\epsilon$ and $\sigma$ affects the potential:

```python
# Compare different epsilon values (depth of well)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

epsilons = [0.5, 1.0, 2.0]
colors_eps = ['lightblue', 'blue', 'darkblue']

for eps, color in zip(epsilons, colors_eps):
    U = lennard_jones(r, eps, sigma_Ar)
    ax1.plot(r, U, color=color, linewidth=2, label=f'ε = {eps}')
    
ax1.axhline(y=0, color='k', linestyle='--', alpha=0.3)
ax1.set_xlim(0.8, 3.5)
ax1.set_ylim(-2.5, 2)
ax1.set_xlabel('Distance r/σ', fontsize=12)
ax1.set_ylabel('Energy U/ε', fontsize=12)
ax1.set_title('Effect of Well Depth (ε)', fontsize=13, fontweight='bold')
ax1.legend(fontsize=11)
ax1.grid(alpha=0.3)

# Compare different sigma values (atom size)
sigmas = [0.8, 1.0, 1.2]
colors_sig = ['coral', 'red', 'darkred']

for sig, color in zip(sigmas, colors_sig):
    U = lennard_jones(r, epsilon_Ar, sig)
    ax2.plot(r, U, color=color, linewidth=2, label=f'σ = {sig}')
    
ax2.axhline(y=0, color='k', linestyle='--', alpha=0.3)
ax2.set_xlim(0.8, 3.5)
ax2.set_ylim(-1.5, 2)
ax2.set_xlabel('Distance r', fontsize=12)
ax2.set_ylabel('Energy U/ε', fontsize=12)
ax2.set_title('Effect of Atom Size (σ)', fontsize=13, fontweight='bold')
ax2.legend(fontsize=11)
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('LJ_parameter_effects.png', dpi=150, bbox_inches='tight')
plt.show()

print("\n" + "="*60)
print("PARAMETER EFFECTS:")
print("="*60)
print("Epsilon (ε) - Controls interaction strength:")
print("  • Larger ε → deeper well → stronger attraction")
print("  • Affects binding energy, not equilibrium distance")
print("  • Noble gases: He (ε≈0.02 kcal/mol) vs Xe (ε≈0.4 kcal/mol)")
print("\nSigma (σ) - Controls atom size:")
print("  • Larger σ → equilibrium at larger distance")
print("  • Shifts entire curve horizontally")
print("  • Noble gases: He (σ≈2.6 Å) vs Xe (σ≈4.1 Å)")
print("="*60)
```

### 2.6 Real Systems: LJ Parameters for Noble Gases

The LJ potential works remarkably well for noble gases. Parameters are typically fitted to reproduce experimental properties like liquid density or vapor pressure.

| Element | σ (Å) | ε (K) | ε (kcal/mol) | ε/k_B (K) |
|---------|-------|-------|--------------|-----------|
| Helium (He) | 2.556 | 10.22 | 0.020 | 10.22 |
| Neon (Ne) | 2.749 | 36.83 | 0.073 | 36.83 |
| Argon (Ar) | 3.405 | 119.8 | 0.238 | 119.8 |
| Krypton (Kr) | 3.624 | 171.0 | 0.340 | 171.0 |
| Xenon (Xe) | 4.063 | 221.0 | 0.439 | 221.0 |

**Note:** Energy units can be expressed in different ways:
- kcal/mol (chemistry convention)
- kJ/mol (SI units)
- K (Kelvin, when divided by Boltzmann constant: ε/k_B)
- Reduced units (normalized to 1)

### 2.7 Limitations of the LJ Potential

While remarkably successful for noble gases, LJ has important limitations:

**Cannot Describe:**
- ❌ Chemical bonds (covalent, ionic)
- ❌ Directional interactions (angles)
- ❌ Charge (electrostatics)
- ❌ Polarization effects
- ❌ Three-body or many-body forces

**Works Best For:**
- ✅ Noble gases (Ar, Ne, Xe)
- ✅ Coarse-grained models (united atom beads)
- ✅ Non-bonded interactions in molecular force fields
- ✅ Simple fluids at low to moderate pressures

**Key Insight:**  
To simulate molecules (not just atoms), we need to add **bonded interactions** on top of the LJ framework. This is what we'll explore next.

---

## 3. Bond Stretching Potentials

### 3.1 The Need for Bonds

Consider a nitrogen molecule (N₂). Two nitrogen atoms are held together by a strong triple bond. If we only used LJ interactions:

1. The atoms would be attracted toward each other (LJ minimum)
2. But they wouldn't maintain a **fixed** bond length
3. At high temperatures, they might even separate completely

We need an additional potential that represents the **covalent bond** connecting the two atoms.

### 3.2 The Harmonic Bond Potential

The simplest model treats a chemical bond like a spring that obeys Hooke's Law:

$$U_{bond}(r) = \frac{1}{2}k_b (r - r_0)^2$$

**Parameters:**
- $k_b$ = bond force constant (energy/length²), typically in kcal/mol/Ų or kJ/mol/Ų
- $r_0$ = equilibrium bond length (Angstroms)

**Physical Interpretation:**
- The bond prefers to stay at length $r_0$
- Stretching OR compressing costs energy (quadratic penalty)
- The stiffer the bond (larger $k_b$), the harder it is to deform

**The Force:**
$$F_{bond}(r) = -\frac{dU_{bond}}{dr} = -k_b(r - r_0)$$

This is the classic spring force: $F = -kx$ where $x = r - r_0$ is the displacement.

### 3.3 Interactive Example: Comparing Bond Strengths

```python
# Bond potential visualization
r_bond = np.linspace(0.5, 2.5, 300)

# Typical bond parameters
r0 = 1.0  # Equilibrium bond length (Angstroms)
k_single = 200   # kcal/mol/Ų - single bond (C-C)
k_double = 500   # kcal/mol/Ų - double bond (C=C)
k_triple = 1000  # kcal/mol/Ų - triple bond (N≡N)

U_single = 0.5 * k_single * (r_bond - r0)**2
U_double = 0.5 * k_double * (r_bond - r0)**2
U_triple = 0.5 * k_triple * (r_bond - r0)**2

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Potential energy
ax1.plot(r_bond, U_single, 'b-', linewidth=2.5, label=f'Single (k={k_single})')
ax1.plot(r_bond, U_double, 'g-', linewidth=2.5, label=f'Double (k={k_double})')
ax1.plot(r_bond, U_triple, 'r-', linewidth=2.5, label=f'Triple (k={k_triple})')
ax1.axvline(x=r0, color='k', linestyle='--', alpha=0.5, label=f'r₀ = {r0} Å')

ax1.set_xlim(0.5, 2.5)
ax1.set_ylim(0, 250)
ax1.set_xlabel('Bond Length r (Å)', fontsize=12)
ax1.set_ylabel('Energy (kcal/mol)', fontsize=12)
ax1.set_title('Harmonic Bond Potential', fontsize=13, fontweight='bold')
ax1.legend(fontsize=10)
ax1.grid(alpha=0.3)

# Force
F_single = -k_single * (r_bond - r0)
F_double = -k_double * (r_bond - r0)
F_triple = -k_triple * (r_bond - r0)

ax2.plot(r_bond, F_single, 'b-', linewidth=2.5, label=f'Single (k={k_single})')
ax2.plot(r_bond, F_double, 'g-', linewidth=2.5, label=f'Double (k={k_double})')
ax2.plot(r_bond, F_triple, 'r-', linewidth=2.5, label=f'Triple (k={k_triple})')
ax2.axhline(y=0, color='k', linestyle='--', alpha=0.3)
ax2.axvline(x=r0, color='k', linestyle='--', alpha=0.5)

ax2.set_xlim(0.5, 2.5)
ax2.set_ylim(-600, 600)
ax2.set_xlabel('Bond Length r (Å)', fontsize=12)
ax2.set_ylabel('Force (kcal/mol/Å)', fontsize=12)
ax2.set_title('Bond Restoring Force', fontsize=13, fontweight='bold')
ax2.legend(fontsize=10)
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('harmonic_bond_potential.png', dpi=150, bbox_inches='tight')
plt.show()

print("\n" + "="*60)
print("BOND CHARACTERISTICS:")
print("="*60)
print(f"Equilibrium length: {r0} Å (no force)")
print(f"\nStretching by 0.1 Å from equilibrium:")
print(f"  Single bond: {0.5 * k_single * 0.1**2:.2f} kcal/mol")
print(f"  Double bond: {0.5 * k_double * 0.1**2:.2f} kcal/mol")
print(f"  Triple bond: {0.5 * k_triple * 0.1**2:.2f} kcal/mol")
print("\n→ Stiffer bonds (larger k) resist deformation more strongly")
print("="*60)
```

### 3.4 Limitations: The Anharmonicity Problem

The harmonic approximation has a critical flaw: **bonds can break in real molecules, but not in the harmonic model**.

Look at the potential energy curve:
- As $r \to 0$: $U \to +\infty$ (infinite compression energy)
- As $r \to \infty$: $U \to +\infty$ (infinite stretching energy)

This means:
- A harmonic bond can never break (no bond dissociation)
- It becomes unphysically stiff at large displacements
- Not suitable for simulating bond-breaking chemistry

**When to Use Harmonic:**
- ✅ Small vibrations around equilibrium (most MD simulations)
- ✅ Systems at moderate temperatures
- ✅ Non-reactive chemistry

**When NOT to Use Harmonic:**
- ❌ Polymer melts (chain entanglement at high strain)
- ❌ Crack propagation in materials
- ❌ Any chemistry involving bond breaking/forming

### 3.5 Advanced Bond Potentials

For systems where bonds can break, more sophisticated potentials exist:

#### FENE Potential (Finitely Extensible Nonlinear Elastic)

Popular in polymer simulations:

$$U_{FENE}(r) = -\frac{1}{2}k_b R_0^2 \ln\left[1 - \left(\frac{r-r_0}{R_0}\right)^2\right]$$

- Has a maximum extension $R_0$ beyond which $U \to \infty$
- Allows large fluctuations without breaking
- Common in coarse-grained polymer models

#### Morse Potential

Realistic potential that allows bond dissociation:

$$U_{Morse}(r) = D_e \left[1 - e^{-\alpha(r-r_0)}\right]^2$$

- $D_e$ = bond dissociation energy
- At $r \to \infty$: $U \to D_e$ (bond breaks, but doesn't diverge)
- Used in studies of chemical reactions

```python
# Compare different bond potentials
r_extended = np.linspace(0.5, 4.0, 500)

# Harmonic (reference)
U_harm = 0.5 * 300 * (r_extended - r0)**2

# FENE parameters
k_fene = 30  # kcal/mol/Ų
R0_fene = 1.5  # Maximum extension
# Avoid log of negative number
r_safe = r_extended[np.abs(r_extended - r0) < R0_fene]
U_fene_safe = -0.5 * k_fene * R0_fene**2 * np.log(1 - ((r_safe - r0)/R0_fene)**2)

# Morse parameters
De = 100  # kcal/mol (dissociation energy)
alpha = 2.0  # Å⁻¹ (width parameter)
U_morse = De * (1 - np.exp(-alpha * (r_extended - r0)))**2

plt.figure(figsize=(10, 6))
plt.plot(r_extended, U_harm, 'b-', linewidth=2.5, label='Harmonic', alpha=0.7)
plt.plot(r_safe, U_fene_safe, 'g-', linewidth=2.5, label='FENE', alpha=0.7)
plt.plot(r_extended, U_morse, 'r-', linewidth=2.5, label='Morse', alpha=0.7)

plt.axvline(x=r0, color='k', linestyle='--', alpha=0.4, label='Equilibrium')
plt.axhline(y=De, color='orange', linestyle=':', linewidth=2, 
            label=f'Dissociation Energy = {De} kcal/mol')

plt.xlim(0.5, 4.0)
plt.ylim(0, 250)
plt.xlabel('Bond Length r (Å)', fontsize=12)
plt.ylabel('Energy (kcal/mol)', fontsize=12)
plt.title('Comparison of Bond Potentials', fontsize=14, fontweight='bold')
plt.legend(fontsize=11)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('bond_potential_comparison.png', dpi=150, bbox_inches='tight')
plt.show()

print("\n" + "="*60)
print("BOND POTENTIAL COMPARISON:")
print("="*60)
print("Harmonic:")
print("  • Diverges at large r → bond cannot break")
print("  • Simple, fast to compute")
print("  • Good for equilibrium simulations")
print("\nFENE:")
print("  • Maximum extension limit")
print("  • Popular in polymer physics")
print("  • Prevents chain crossing")
print("\nMorse:")
print("  • Realistic dissociation behavior")
print("  • More expensive to compute")
print("  • Used in reactive simulations")
print("="*60)
```

---

## 4. Angle Bending Potentials

### 4.1 Why Molecules Have Shapes

Atoms in molecules don't just bond—they bond at **specific angles**. Consider water (H₂O):
- The two O-H bonds have a specific angle: 104.5°
- This angle is not arbitrary—it arises from the sp³ hybridization of oxygen's electron orbitals

Without angle potentials, a water molecule in MD would be floppy, with the H-O-H angle fluctuating wildly. We need a potential that maintains molecular geometry.

### 4.2 The Harmonic Angle Potential

Similar to bonds, we model angle bending as a spring:

$$U_{angle}(\theta) = \frac{1}{2}k_\theta (\theta - \theta_0)^2$$

**Parameters:**
- $k_\theta$ = angle force constant (energy/radian²), typically in kcal/mol/rad²
- $\theta_0$ = equilibrium angle (degrees or radians)

**Key Points:**
- Three atoms are involved: i-j-k (j is the central atom)
- The angle is measured between bonds j-i and j-k
- Energy is minimized when $\theta = \theta_0$

**The Force (Torque):**
The derivative with respect to angle gives a torque:
$$\tau = -\frac{dU_{angle}}{d\theta} = -k_\theta(\theta - \theta_0)$$

### 4.3 Common Molecular Geometries

Different hybridizations lead to different preferred angles:

| Geometry | Hybridization | Ideal Angle | Example |
|----------|---------------|-------------|---------|
| Linear | sp | 180° | CO₂, HCN |
| Trigonal planar | sp² | 120° | Benzene, ethene |
| Tetrahedral | sp³ | 109.47° | Methane (CH₄), diamond |
| Bent (water) | sp³ | 104.5° | H₂O |
| Trigonal pyramidal | sp³ | ~107° | NH₃ (ammonia) |

```python
# Visualize angle potential for different molecular geometries
theta_deg = np.linspace(60, 180, 300)
theta_rad = np.deg2rad(theta_deg)

# Different angle force constants (typical values)
k_soft = 50    # kcal/mol/rad² - flexible angle
k_medium = 100 # kcal/mol/rad² - standard angle
k_stiff = 200  # kcal/mol/rad² - rigid angle

# Different equilibrium angles
theta0_linear = np.deg2rad(180)  # CO₂
theta0_trigonal = np.deg2rad(120)  # Benzene
theta0_tetrahedral = np.deg2rad(109.47)  # Methane
theta0_bent = np.deg2rad(104.5)  # Water

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# --- Plot 1: Effect of force constant ---
U_soft = 0.5 * k_soft * (theta_rad - theta0_tetrahedral)**2
U_medium = 0.5 * k_medium * (theta_rad - theta0_tetrahedral)**2
U_stiff = 0.5 * k_stiff * (theta_rad - theta0_tetrahedral)**2

ax1.plot(theta_deg, U_soft, 'lightblue', linewidth=2.5, label=f'k = {k_soft} (flexible)')
ax1.plot(theta_deg, U_medium, 'blue', linewidth=2.5, label=f'k = {k_medium} (standard)')
ax1.plot(theta_deg, U_stiff, 'darkblue', linewidth=2.5, label=f'k = {k_stiff} (stiff)')
ax1.axvline(x=109.47, color='k', linestyle='--', alpha=0.5, label='θ₀ = 109.47° (sp³)')

ax1.set_xlabel('Angle θ (degrees)', fontsize=12)
ax1.set_ylabel('Energy (kcal/mol)', fontsize=12)
ax1.set_title('Effect of Angle Force Constant', fontsize=13, fontweight='bold')
ax1.set_ylim(0, 15)
ax1.legend(fontsize=10)
ax1.grid(alpha=0.3)

# --- Plot 2: Different equilibrium angles ---
U_linear = 0.5 * k_medium * (theta_rad - theta0_linear)**2
U_trigonal = 0.5 * k_medium * (theta_rad - theta0_trigonal)**2
U_tetrahedral = 0.5 * k_medium * (theta_rad - theta0_tetrahedral)**2
U_bent = 0.5 * k_medium * (theta_rad - theta0_bent)**2

ax2.plot(theta_deg, U_linear, 'purple', linewidth=2.5, label='Linear (180°, CO₂)')
ax2.plot(theta_deg, U_trigonal, 'green', linewidth=2.5, label='Trigonal (120°, benzene)')
ax2.plot(theta_deg, U_tetrahedral, 'blue', linewidth=2.5, label='Tetrahedral (109.5°, CH₄)')
ax2.plot(theta_deg, U_bent, 'red', linewidth=2.5, label='Bent (104.5°, H₂O)')

ax2.set_xlabel('Angle θ (degrees)', fontsize=12)
ax2.set_ylabel('Energy (kcal/mol)', fontsize=12)
ax2.set_title('Different Molecular Geometries', fontsize=13, fontweight='bold')
ax2.set_ylim(0, 15)
ax2.legend(fontsize=10)
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('angle_potentials.png', dpi=150, bbox_inches='tight')
plt.show()

print("\n" + "="*60)
print("ANGLE BENDING EXAMPLES:")
print("="*60)
print("Water (H₂O):")
print(f"  θ₀ = 104.5°, k_θ ≈ 100 kcal/mol/rad²")
print(f"  Bending by 10° costs: {0.5 * 100 * np.deg2rad(10)**2:.2f} kcal/mol")
print("\nMethane (CH₄):")
print(f"  θ₀ = 109.47° (perfect tetrahedral)")
print(f"  All H-C-H angles equal in equilibrium")
print("="*60)
```

### 4.4 When Angles Matter: A Physical Example

**Case Study: Why Ice Floats**

Water's unusual property (solid less dense than liquid) arises from its angle potential:

1. The 104.5° H-O-H angle forces water molecules into a specific tetrahedral arrangement
2. When water freezes, this geometry creates an **open hexagonal lattice** with empty space
3. Liquid water has more collapsed, disordered angles → higher density
4. Therefore: Ice density (0.92 g/cm³) < Water density (1.00 g/cm³)

If water had a flexible angle (small $k_\theta$), ice would be denser and would sink—fundamentally changing Earth's climate and the evolution of life.

---


## 5. Dihedral (Torsional) Potentials

### 5.1 Rotation Around Bonds: The Fourth Dimension of Molecular Structure

So far we've controlled:
- **Bond length** (distance between two atoms)
- **Angle** (geometry between three atoms)

But molecules have one more degree of freedom: **rotation around bonds**. Consider ethane (C₂H₆):
- The two methyl groups (CH₃) can rotate relative to each other
- Different rotational configurations have different energies
- This rotation is described by the **dihedral angle** (also called torsion angle)

### 5.2 Defining the Dihedral Angle

A dihedral angle involves **four atoms in sequence**: i-j-k-l

**Geometric Definition:**
- Form two planes:
  - Plane 1: atoms i, j, k
  - Plane 2: atoms j, k, l
- The dihedral angle φ is the angle between these two planes

**Visualization:**
- φ = 0° : all four atoms are coplanar (cis configuration)
- φ = 180° : atoms i and l are on opposite sides (trans configuration)
- φ = ±90° : perpendicular planes

```python
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Simple visualization of dihedral angle concept
fig = plt.figure(figsize=(14, 5))

# Three different dihedral angles
angles = [0, 90, 180]
titles = ['φ = 0° (Cis/Eclipsed)', 'φ = 90° (Gauche)', 'φ = 180° (Trans/Staggered)']

for idx, (phi, title) in enumerate(zip(angles, titles)):
    ax = fig.add_subplot(1, 3, idx+1, projection='3d')
    
    # Define four atoms in a chain (simplified ethane-like)
    # Atom positions for visualization
    atom1 = np.array([0, 0, 0])
    atom2 = np.array([1, 0, 0])
    atom3 = np.array([2, 0, 0])
    # Atom 4 rotates around the j-k bond
    phi_rad = np.deg2rad(phi)
    atom4 = np.array([3, np.cos(phi_rad), np.sin(phi_rad)])
    
    atoms = np.array([atom1, atom2, atom3, atom4])
    
    # Plot atoms
    ax.scatter(atoms[:, 0], atoms[:, 1], atoms[:, 2], 
               c=['red', 'blue', 'blue', 'green'], s=200, alpha=0.8)
    
    # Plot bonds
    for i in range(3):
        ax.plot([atoms[i, 0], atoms[i+1, 0]], 
                [atoms[i, 1], atoms[i+1, 1]], 
                [atoms[i, 2], atoms[i+1, 2]], 'k-', linewidth=2)
    
    # Labels
    ax.text(atom1[0], atom1[1], atom1[2]-0.3, 'i', fontsize=12, ha='center')
    ax.text(atom2[0], atom2[1], atom2[2]-0.3, 'j', fontsize=12, ha='center')
    ax.text(atom3[0], atom3[1], atom3[2]-0.3, 'k', fontsize=12, ha='center')
    ax.text(atom4[0], atom4[1], atom4[2]-0.3, 'l', fontsize=12, ha='center')
    
    ax.set_xlim(-0.5, 3.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_zlim(-1.5, 1.5)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title, fontsize=11, fontweight='bold')

plt.tight_layout()
plt.savefig('dihedral_angle_definition.png', dpi=150, bbox_inches='tight')
plt.show()

print("\n" + "="*60)
print("DIHEDRAL ANGLE CONVENTION:")
print("="*60)
print("• Four atoms: i-j-k-l (bonded in sequence)")
print("• Central bond: j-k (axis of rotation)")
print("• φ = 0°: cis/eclipsed (i and l on same side)")
print("• φ = 180°: trans/staggered (i and l on opposite sides)")
print("• Positive rotation: right-hand rule around j→k vector")
print("="*60)
```

### 5.3 The Cosine Series Potential

Unlike bonds and angles (harmonic), dihedrals are **periodic**—rotating by 360° brings you back to the same configuration. The standard form is a Fourier series:

$$U_{dihedral}(\phi) = \sum_{n=1}^{N} \frac{V_n}{2} \left[1 + \cos(n\phi - \phi_n)\right]$$

**Parameters:**
- $V_n$ = barrier height for the n-th term (kcal/mol)
- $n$ = multiplicity (number of minima in 360°)
- $\phi_n$ = phase shift (determines where minima occur)

**Common Special Cases:**

**1-fold (n=1):** One minimum per 360° rotation
$$U_1(\phi) = \frac{V_1}{2}[1 + \cos(\phi)]$$
- Used for maintaining planarity (e.g., peptide bonds)

**2-fold (n=2):** Two minima per 360° rotation  
$$U_2(\phi) = \frac{V_2}{2}[1 + \cos(2\phi)]$$

**3-fold (n=3):** Three minima per 360° rotation
$$U_3(\phi) = \frac{V_3}{2}[1 + \cos(3\phi)]$$
- Classic for C-C bond rotation in hydrocarbons

### 5.4 Interactive Example: Ethane Rotation

Ethane (H₃C-CH₃) is the prototypical example. Rotation around the C-C bond shows a **3-fold barrier**:

```python
# Dihedral potential for ethane rotation
phi_deg = np.linspace(0, 360, 500)
phi_rad = np.deg2rad(phi_deg)

# Ethane C-C rotation parameters (typical values)
V1 = 0.0      # No 1-fold component
V2 = 0.0      # No 2-fold component
V3 = 1.4      # 3-fold barrier height (kcal/mol)

# Calculate potential
U_ethane = (V3/2) * (1 + np.cos(3*phi_rad))

# Also show individual n-fold components for comparison
U_1fold = (2.0/2) * (1 + np.cos(phi_rad))      # Example 1-fold
U_2fold = (1.5/2) * (1 + np.cos(2*phi_rad))    # Example 2-fold
U_3fold = (V3/2) * (1 + np.cos(3*phi_rad))     # Ethane

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# --- Plot 1: n-fold comparison ---
ax1.plot(phi_deg, U_1fold, 'b-', linewidth=2.5, label='n=1 (1 minimum)')
ax1.plot(phi_deg, U_2fold, 'g-', linewidth=2.5, label='n=2 (2 minima)')
ax1.plot(phi_deg, U_3fold, 'r-', linewidth=2.5, label='n=3 (3 minima)')

ax1.set_xlabel('Dihedral Angle φ (degrees)', fontsize=12)
ax1.set_ylabel('Energy (kcal/mol)', fontsize=12)
ax1.set_title('Effect of Multiplicity (n)', fontsize=13, fontweight='bold')
ax1.set_xlim(0, 360)
ax1.set_xticks([0, 60, 120, 180, 240, 300, 360])
ax1.legend(fontsize=11)
ax1.grid(alpha=0.3)

# --- Plot 2: Ethane with annotations ---
ax2.plot(phi_deg, U_ethane, 'darkblue', linewidth=3, label='Ethane C-C rotation')

# Mark special conformations
eclipsed_angles = [0, 120, 240, 360]
staggered_angles = [60, 180, 300]

for angle in eclipsed_angles:
    U_val = (V3/2) * (1 + np.cos(3*np.deg2rad(angle)))
    ax2.plot(angle, U_val, 'ro', markersize=10)
    if angle == 0:
        ax2.annotate('Eclipsed\n(max E)', xy=(angle, U_val), 
                    xytext=(angle, U_val+0.3), fontsize=9, ha='center',
                    arrowprops=dict(arrowstyle='->', color='red', lw=1.5))

for angle in staggered_angles:
    U_val = (V3/2) * (1 + np.cos(3*np.deg2rad(angle)))
    ax2.plot(angle, U_val, 'go', markersize=10)
    if angle == 60:
        ax2.annotate('Staggered\n(min E)', xy=(angle, U_val), 
                    xytext=(angle, U_val-0.4), fontsize=9, ha='center',
                    arrowprops=dict(arrowstyle='->', color='green', lw=1.5))

ax2.axhline(y=0, color='k', linestyle='--', alpha=0.3)
ax2.set_xlabel('Dihedral Angle φ (degrees)', fontsize=12)
ax2.set_ylabel('Energy (kcal/mol)', fontsize=12)
ax2.set_title('Ethane: 3-fold Rotational Barrier', fontsize=13, fontweight='bold')
ax2.set_xlim(0, 360)
ax2.set_ylim(-0.2, 1.6)
ax2.set_xticks([0, 60, 120, 180, 240, 300, 360])
ax2.legend(fontsize=11)
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('dihedral_ethane.png', dpi=150, bbox_inches='tight')
plt.show()

print("\n" + "="*60)
print("ETHANE ROTATIONAL BARRIER:")
print("="*60)
print(f"Barrier height: {V3} kcal/mol")
print("\nStaggered conformations (φ = 60°, 180°, 300°):")
print(f"  Energy: {(V3/2) * (1 + np.cos(3*np.deg2rad(60))):.3f} kcal/mol (minimum)")
print("\nEclipsed conformations (φ = 0°, 120°, 240°):")
print(f"  Energy: {(V3/2) * (1 + np.cos(3*np.deg2rad(0))):.3f} kcal/mol (maximum)")
print("\nPhysical origin:")
print("  • Staggered: H atoms on opposite carbons are far apart")
print("  • Eclipsed: H atoms align, causing steric repulsion")
print("="*60)
```

### 5.5 Combining Multiple Terms

Real molecules often need multiple n-fold terms. For example, peptide bonds (C-N in proteins) use:

$$U_{peptide}(\phi) = V_1[1 + \cos(\phi)] + V_2[1 + \cos(2\phi)]$$

This creates an asymmetric potential that strongly favors planarity (trans configuration at φ ≈ 180°).

```python
# Example: Peptide bond with multiple terms
V1_peptide = 3.0   # Strong 1-fold (maintains planarity)
V2_peptide = 1.0   # Weak 2-fold (fine-tuning)

U_peptide = (V1_peptide/2) * (1 + np.cos(phi_rad)) + \
            (V2_peptide/2) * (1 + np.cos(2*phi_rad))

plt.figure(figsize=(10, 6))
plt.plot(phi_deg, U_peptide, 'purple', linewidth=3, label='Peptide bond (n=1+2)')
plt.axvline(x=180, color='green', linestyle='--', linewidth=2, 
            alpha=0.7, label='Trans (φ=180°, preferred)')
plt.axvline(x=0, color='red', linestyle='--', linewidth=2, 
            alpha=0.7, label='Cis (φ=0°, rare)')

plt.xlabel('Dihedral Angle φ (degrees)', fontsize=12)
plt.ylabel('Energy (kcal/mol)', fontsize=12)
plt.title('Peptide Bond: Asymmetric Barrier (Combined n=1+2)', fontsize=13, fontweight='bold')
plt.xlim(0, 360)
plt.xticks([0, 90, 180, 270, 360])
plt.legend(fontsize=11)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('dihedral_peptide.png', dpi=150, bbox_inches='tight')
plt.show()

print("\n" + "="*60)
print("PEPTIDE BOND CHARACTERISTICS:")
print("="*60)
print("Trans (φ ≈ 180°):")
print(f"  Energy: {(V1_peptide/2)*(1+np.cos(np.deg2rad(180))) + (V2_peptide/2)*(1+np.cos(2*np.deg2rad(180))):.2f} kcal/mol")
print(f"  Stability: Strongly preferred (>99% in proteins)")
print("\nCis (φ ≈ 0°):")
print(f"  Energy: {(V1_peptide/2)*(1+np.cos(np.deg2rad(0))) + (V2_peptide/2)*(1+np.cos(2*np.deg2rad(0))):.2f} kcal/mol")
print(f"  Stability: Rare except before proline residues")
print("\nBarrier (~7 kcal/mol) prevents spontaneous isomerization")
print("="*60)
```

### 5.6 Why Dihedrals Matter: Protein Folding

The **Ramachandran plot** in protein biochemistry is fundamentally a map of dihedral angle energy:
- φ (phi): rotation around N-Cα bond
- ψ (psi): rotation around Cα-C bond

Allowed regions of the plot correspond to low-energy dihedral combinations. Without proper dihedral potentials, proteins would unfold completely in MD simulations.

---

## 6. Electrostatic Interactions

### 6.1 The Long-Range Force

All the potentials we've discussed so far are **short-range**:
- LJ: significant only within ~10 Å
- Bonds/Angles/Dihedrals: operate on directly connected atoms

But many molecules have **charged atoms** or **partial charges** (polar bonds). The electrostatic interaction between charges is **long-range** and follows Coulomb's law:

$$U_{elec}(r) = \frac{1}{4\pi\epsilon_0} \frac{q_i q_j}{r_{ij}} = k_e \frac{q_i q_j}{r_{ij}}$$

**Parameters:**
- $q_i, q_j$ = partial charges on atoms i and j (in units of electron charge, e)
- $r_{ij}$ = distance between atoms
- $k_e = \frac{1}{4\pi\epsilon_0} = 332.06$ kcal·Å/(mol·e²) in common MD units
- $\epsilon_0$ = permittivity of free space

**Key Difference from LJ:**
- LJ decays as $r^{-6}$ → negligible beyond cutoff (~10 Å)
- Coulomb decays as $r^{-1}$ → significant even at 50-100 Å
- This makes electrostatics computationally expensive

### 6.2 Partial Charges: The Concept

In a covalent bond between different elements, electrons are not shared equally:
- More electronegative atom (e.g., O, N, F) pulls electrons toward itself
- Becomes partially negative (δ-)
- Partner atom becomes partially positive (δ+)

**Example: Water (H₂O)**
- Oxygen is highly electronegative
- Each O-H bond is polar
- Typical partial charges:
  - O: q = -0.82 e
  - H: q = +0.41 e (each)
  - Total: -0.82 + 2(0.41) = 0 (molecule is neutral)

### 6.3 Interactive Example: Charge-Charge Interactions

```python
# Coulomb potential visualization
r_elec = np.linspace(1, 20, 300)  # Distance in Angstroms

# Coulomb constant in MD units (kcal·Å/mol/e²)
k_e = 332.06

# Different charge combinations
# Like charges (repulsion)
q_plus_plus = k_e * (+1) * (+1) / r_elec
q_minus_minus = k_e * (-1) * (-1) / r_elec

# Opposite charges (attraction)  
q_plus_minus = k_e * (+1) * (-1) / r_elec

# Partial charges (water H-H)
q_partial = k_e * (+0.41) * (+0.41) / r_elec

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# --- Plot 1: Full charges ---
ax1.plot(r_elec, q_plus_plus, 'r-', linewidth=2.5, label='(+1)(+1) repulsion')
ax1.plot(r_elec, q_minus_minus, 'orange', linewidth=2.5, label='(-1)(-1) repulsion')
ax1.plot(r_elec, q_plus_minus, 'b-', linewidth=2.5, label='(+1)(-1) attraction')
ax1.axhline(y=0, color='k', linestyle='--', alpha=0.3)

ax1.set_xlabel('Distance r (Å)', fontsize=12)
ax1.set_ylabel('Coulomb Energy (kcal/mol)', fontsize=12)
ax1.set_title('Full Charges (±1 e)', fontsize=13, fontweight='bold')
ax1.set_xlim(1, 20)
ax1.set_ylim(-350, 350)
ax1.legend(fontsize=11)
ax1.grid(alpha=0.3)

# --- Plot 2: Comparison with LJ ---
# Use Argon LJ parameters for comparison
epsilon_Ar = 0.238  # kcal/mol
sigma_Ar = 3.405    # Angstrom
U_LJ = 4 * epsilon_Ar * ((sigma_Ar/r_elec)**12 - (sigma_Ar/r_elec)**6)

ax2.plot(r_elec, q_partial, 'purple', linewidth=2.5, label='Coulomb: (+0.41)(+0.41)')
ax2.plot(r_elec, U_LJ, 'green', linewidth=2.5, label='LJ: Argon', linestyle='--')
ax2.axhline(y=0, color='k', linestyle='--', alpha=0.3)

ax2.set_xlabel('Distance r (Å)', fontsize=12)
ax2.set_ylabel('Energy (kcal/mol)', fontsize=12)
ax2.set_title('Electrostatic vs. van der Waals', fontsize=13, fontweight='bold')
ax2.set_xlim(1, 20)
ax2.set_ylim(-1, 3)
ax2.legend(fontsize=11)
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('coulomb_interactions.png', dpi=150, bbox_inches='tight')
plt.show()

print("\n" + "="*60)
print("ELECTROSTATIC ENERGY SCALE:")
print("="*60)
print(f"Full charges (+1, -1) at 5 Å:")
print(f"  U = {k_e * 1 * -1 / 5:.1f} kcal/mol (VERY strong)")
print(f"\nPartial charges (H₂O: +0.41, -0.82) at 3 Å:")
print(f"  U = {k_e * 0.41 * -0.82 / 3:.1f} kcal/mol (moderate)")
print(f"\nCompare to:")
print(f"  • Thermal energy (RT): ~0.6 kcal/mol at 300K")
print(f"  • LJ well depth (Ar): 0.24 kcal/mol")
print(f"  • Hydrogen bond: ~5 kcal/mol")
print("\n→ Electrostatics dominate charged/polar systems")
print("="*60)
```

### 6.4 Combining LJ and Coulomb: The Non-Bonded Potential

In molecular force fields, atoms interact through BOTH LJ and Coulomb:

$$U_{non-bonded}(r_{ij}) = 4\epsilon_{ij}\left[\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{r_{ij}}\right)^6\right] + k_e\frac{q_i q_j}{r_{ij}}$$

**Combining Rules for LJ:**
When two different atom types interact (e.g., O-H rather than O-O), we need to combine their individual parameters:

**Lorentz-Berthelot (most common):**
- $\sigma_{ij} = \frac{\sigma_i + \sigma_j}{2}$ (arithmetic mean)
- $\epsilon_{ij} = \sqrt{\epsilon_i \epsilon_j}$ (geometric mean)

**Geometric mean for both (alternative):**
- $\sigma_{ij} = \sqrt{\sigma_i \sigma_j}$
- $\epsilon_{ij} = \sqrt{\epsilon_i \epsilon_j}$

```python
# Example: O-H interaction in water
# Oxygen parameters (TIP3P water model)
sigma_O = 3.1507   # Angstrom
epsilon_O = 0.1521  # kcal/mol
q_O = -0.834       # electron charge

# Hydrogen parameters (TIP3P)
sigma_H = 0.0      # Hydrogen is massless point charge in TIP3P
epsilon_H = 0.0
q_H = 0.417

# Combined parameters (Lorentz-Berthelot)
sigma_OH = (sigma_O + sigma_H) / 2
epsilon_OH = np.sqrt(epsilon_O * epsilon_H)

r_water = np.linspace(0.8, 8, 300)

# Calculate O-H interaction
U_LJ_OH = 4 * epsilon_OH * ((sigma_OH/r_water)**12 - (sigma_OH/r_water)**6)
U_Coul_OH = k_e * q_O * q_H / r_water
U_total_OH = U_LJ_OH + U_Coul_OH

# Calculate O-O interaction for comparison
sigma_OO = sigma_O
epsilon_OO = epsilon_O
U_LJ_OO = 4 * epsilon_OO * ((sigma_OO/r_water)**12 - (sigma_OO/r_water)**6)
U_Coul_OO = k_e * q_O * q_O / r_water
U_total_OO = U_LJ_OO + U_Coul_OO

plt.figure(figsize=(10, 6))

plt.plot(r_water, U_LJ_OO, 'b--', linewidth=2, alpha=0.6, label='O-O: LJ only')
plt.plot(r_water, U_Coul_OO, 'b:', linewidth=2, alpha=0.6, label='O-O: Coulomb only')
plt.plot(r_water, U_total_OO, 'b-', linewidth=2.5, label='O-O: Total')

plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
plt.xlabel('Distance r (Å)', fontsize=12)
plt.ylabel('Energy (kcal/mol)', fontsize=12)
plt.title('Water: Non-Bonded O-O Interaction (TIP3P model)', fontsize=13, fontweight='bold')
plt.xlim(0.8, 8)
plt.ylim(-10, 20)
plt.legend(fontsize=11)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('water_nonbonded.png', dpi=150, bbox_inches='tight')
plt.show()

print("\n" + "="*60)
print("WATER OXYGEN-OXYGEN INTERACTION:")
print("="*60)
print(f"At typical O-O distance (2.8 Å):")
print(f"  LJ contribution: {4*epsilon_OO*((sigma_OO/2.8)**12-(sigma_OO/2.8)**6):.3f} kcal/mol")
print(f"  Coulomb contribution: {k_e*q_O*q_O/2.8:.3f} kcal/mol")
print(f"  Total: {4*epsilon_OO*((sigma_OO/2.8)**12-(sigma_OO/2.8)**6) + k_e*q_O*q_O/2.8:.3f} kcal/mol")
print("\n→ Electrostatics provide strong repulsion (same sign)")
print("→ LJ provides weak attraction at equilibrium")
print("="*60)
```

### 6.5 The Computational Challenge: Long-Range Electrostatics

Because Coulomb interactions decay slowly ($r^{-1}$), we cannot use a simple cutoff like we do for LJ. Special algorithms are required:

**Methods:**
1. **Ewald Summation** - Exact, but slow ($O(N^{3/2})$)
2. **Particle Mesh Ewald (PME)** - Fast approximation ($O(N \log N)$)
3. **Reaction Field** - Continuum approximation beyond cutoff

We'll discuss these in detail when we cover LAMMPS implementation in the tutorials.

---

## 7. The Complete Force Field: Putting It All Together

### 7.1 Total Potential Energy Function

A modern molecular force field combines all the terms we've discussed:

$$U_{total} = U_{bonded} + U_{non-bonded}$$

**Bonded Terms** (intramolecular):
$$U_{bonded} = \sum_{bonds} U_{bond} + \sum_{angles} U_{angle} + \sum_{dihedrals} U_{dihedral} + \sum_{impropers} U_{improper}$$

**Non-Bonded Terms** (inter + intramolecular):
$$U_{non-bonded} = \sum_{i<j} \left[ U_{LJ}(r_{ij}) + U_{Coulomb}(r_{ij}) \right]$$

### 7.2 Force Calculation: The Derivative

Once we have $U_{total}(\mathbf{r}^N)$, the force on each atom is:

$$\mathbf{F}_i = -\nabla_i U_{total} = -\frac{\partial U_{total}}{\partial \mathbf{r}_i}$$

Each component contributes to the force:
- **Bond force**: Acts along bond direction, magnitude $\propto (r - r_0)$
- **Angle force**: Acts perpendicular to bond, creating torque
- **Dihedral force**: Complex, affects all four atoms
- **LJ force**: Radial, from pairwise distances
- **Coulomb force**: Radial, from all charged pairs

### 7.3 Exclusions and Special Rules

**Important:** Not all atom pairs use the non-bonded potential.

**Standard Exclusion Rules:**
1. **1-2 pairs** (directly bonded): Non-bonded **excluded** (bond potential used instead)
2. **1-3 pairs** (separated by 2 bonds): Non-bonded **excluded** (angle potential connects them)
3. **1-4 pairs** (separated by 3 bonds): Non-bonded **scaled** (typically 50% strength)
4. **1-5+ pairs**: Full non-bonded interaction

**Why?**
- Double-counting: The bond potential already describes 1-2 interaction
- Overpolarization: Full LJ+Coulomb on 1-3 would be too strong
- 1-4 scaling: Empirical fitting to reproduce conformational energies

```python
# Visualization of exclusion rules in a simple molecule
fig, ax = plt.subplots(figsize=(10, 6))

# Draw a simple hydrocarbon chain
positions = {
    'C1': (1, 2),
    'C2': (3, 2),
    'C3': (5, 2),
    'C4': (7, 2),
    'H1': (1, 3.5),
    'H2': (3, 3.5),
    'H3': (5, 3.5),
}

# Draw atoms
for name, pos in positions.items():
    color = 'gray' if 'C' in name else 'lightblue'
    ax.scatter(*pos, s=800, c=color, edgecolors='black', linewidths=2, zorder=3)
    ax.text(pos[0], pos[1], name, ha='center', va='center', fontsize=11, fontweight='bold')

# Draw bonds
bonds = [('C1', 'C2'), ('C2', 'C3'), ('C3', 'C4'), 
         ('C1', 'H1'), ('C2', 'H2'), ('C3', 'H3')]
for atom1, atom2 in bonds:
    pos1, pos2 = positions[atom1], positions[atom2]
    ax.plot([pos1[0], pos2[0]], [pos1[1], pos2[1]], 'k-', linewidth=3, zorder=1)

# Annotate different pair types
# 1-2 (bonded)
ax.annotate('', xy=positions['C2'], xytext=positions['C1'],
            arrowprops=dict(arrowstyle='<->', color='red', lw=2.5))
ax.text(2, 1.3, '1-2: Excluded', ha='center', fontsize=10, color='red', fontweight='bold')

# 1-3
ax.annotate('', xy=positions['C3'], xytext=positions['C1'],
            arrowprops=dict(arrowstyle='<->', color='orange', lw=2, linestyle='--'))
ax.text(3, 0.8, '1-3: Excluded', ha='center', fontsize=10, color='orange', fontweight='bold')

# 1-4
ax.annotate('', xy=positions['C4'], xytext=positions['C1'],
            arrowprops=dict(arrowstyle='<->', color='blue', lw=2, linestyle=':'))
ax.text(4, 0.3, '1-4: Scaled (50%)', ha='center', fontsize=10, color='blue', fontweight='bold')

ax.set_xlim(0, 8)
ax.set_ylim(0, 4)
ax.axis('off')
ax.set_title('Non-Bonded Exclusion Rules', fontsize=14, fontweight='bold', pad=20)

plt.tight_layout()
plt.savefig('exclusion_rules.png', dpi=150, bbox_inches='tight')
plt.show()

print("\n" + "="*60)
print("EXCLUSION RULES SUMMARY:")
print("="*60)
print("Pair Type  | Separation | LJ + Coulomb")
print("-----------|------------|-------------")
print("1-2        | 1 bond     | Excluded (use bond potential)")
print("1-3        | 2 bonds    | Excluded (use angle potential)")
print("1-4        | 3 bonds    | Scaled (typically 0.5×)")
print("1-5+       | 4+ bonds   | Full interaction")
print("="*60)
```

### 7.4 Computational Cost Breakdown

For a typical protein simulation with N atoms:

| Component | Scaling | % of Time | Notes |
|-----------|---------|-----------|-------|
| **Bond stretching** | O(N) | ~2% | Few bonds per atom |
| **Angle bending** | O(N) | ~3% | Few angles per atom |
| **Dihedrals** | O(N) | ~5% | Few dihedrals per atom |
| **LJ (cutoff)** | O(N) | ~15% | With neighbor lists |
| **Coulomb (PME)** | O(N log N) | ~70% | Long-range, expensive |
| **Integration** | O(N) | ~5% | Velocity Verlet |

**Key Insight:**  
Electrostatics dominate the computational cost. Efficient PME implementation is critical for performance.

### 7.5 Summary: The Hierarchy Revisited

We've built up from simple to complex:

```
Level 1: Noble gases
  └─ LJ only

Level 2: Diatomic molecules (N₂, O₂)
  └─ LJ + bonds

Level 3: Small molecules (H₂O, CO₂)
  └─ LJ + bonds + angles + charges

Level 4: Organic molecules (ethane, benzene)
  └─ LJ + bonds + angles + dihedrals + charges

Level 5: Proteins, polymers
  └─ Full force field (all terms)
```

Each level captures additional physics, but at the cost of more parameters and complexity.

---

## 8. Practical Considerations

### 8.1 Parameter Sources

Where do these parameters come from?

**Standard Force Fields:**
- **AMBER** (Assisted Model Building with Energy Refinement) - proteins, nucleic acids
- **CHARMM** (Chemistry at Harvard Macromolecular Mechanics) - biomolecules
- **OPLS** (Optimized Potentials for Liquid Simulations) - organic liquids
- **GROMOS** (Groningen Molecular Simulation) - biomolecules
- **TraPPE** (Transferable Potentials for Phase Equilibria) - small molecules

Each force field is a **carefully fitted** set of parameters designed to reproduce specific experimental properties.

**Fitting Strategy:**
1. Choose target properties (e.g., liquid density, heat of vaporization, bond lengths)
2. Run simulations with trial parameters
3. Adjust parameters to minimize error
4. Validate against additional experiments

**Critical:** Don't mix parameters from different force fields—they're fitted as internally consistent sets.

### 8.2 Transferability

A "good" force field has **transferable** parameters:
- Same C-C bond in methane, ethane, and benzene
- Same O-H bond in water, methanol, and sugars

This allows you to simulate new molecules by combining existing atom types.

**Limitations:**
- Parameters fitted for one environment may fail in another
- AMBER protein parameters optimized for aqueous solution may not work in membrane
- Small molecule parameters from gas phase may not describe liquid accurately

### 8.3 Validation: Does Your Force Field Work?

Before trusting MD results, always validate:

**Structural Properties:**
- Bond lengths, angles (compare to X-ray crystallography)
- Radial distribution functions (compare to neutron scattering)
- Density (compare to experiment)

**Thermodynamic Properties:**
- Heat of vaporization
- Free energy of solvation
- Heat capacity

**Dynamic Properties:**
- Diffusion coefficients
- Viscosity
- Relaxation times

If your force field fails validation, your MD results are meaningless.

---

## 9. Next Steps: Advanced Potentials

The classical potentials covered here work well for **non-reactive** systems where bonds don't break or form. In the next lecture, we'll explore:

1. **Complete Force Field Packages** (AMBER, CHARMM, OPLS)
   - How to assign atom types
   - Parameter files and topology generation
   
2. **Many-Body Potentials** (EAM for metals)
   - Why pair potentials fail for metals
   - Embedded Atom Method formalism

3. **Reactive Potentials** (ReaxFF)
   - Bond-order concept
   - Charge equilibration
   - Simulating combustion and catalysis

---

## 10. Key Takeaways

**What You Should Remember:**

1. **MD potentials approximate quantum mechanics** using simple functional forms

2. **The hierarchy of complexity:**
   - LJ: van der Waals forces (noble gases)
   - Bonds: harmonic springs (molecules)
   - Angles: maintain geometry
   - Dihedrals: rotational barriers
   - Coulomb: electrostatics (polar/charged)

3. **Each term has a clear physical meaning:**
   - Force = $-\frac{dU}{dr}$ (always downhill on energy surface)
   - Parameters fitted to reproduce experiments

4. **Computational cost matters:**
   - Electrostatics dominate (70% of CPU time)
   - Smart algorithms (PME) essential

5. **Force field validation is critical:**
   - Never trust results without checking against experiments
   - Transferability is limited—use appropriate force field for your system

**In the tutorials, you'll see:**
- How to implement these potentials in LAMMPS
- Progressive examples showing what each term adds
- How choice of potential affects predicted properties