# Energetics and Free Energy

---

## Part 1: Static Energetics

### 1. The Potential Energy Surface as a Map

Throughout this course, we have used the potential energy surface (PES) as a dynamic object — the landscape that determines forces and accelerates atoms during an MD simulation. But the PES is also a rich static object that can be characterized directly, without running dynamics at all.

Think of the PES as a mountainous terrain in a very high-dimensional space (one dimension per atomic coordinate). MD is like a ball rolling and bouncing across this terrain driven by thermal energy. But we can also ask quieter questions:

- Where are the valleys? (Stable structures)
- How deep are they? (Stability and cohesive energy)
- Which valley is lowest? (Phase stability)
- How wide are they? (Mechanical stiffness)

Answering these questions is the subject of **static energetics** — using a potential to compute energy differences that correspond to real, measurable material properties.

---

### 2. Energy Minimization

The first task is always to find a valley bottom: a configuration where every atom sits at a local minimum of the PES, meaning no net force acts on any atom.

$$\mathbf{F}_i = -\nabla_i U = 0 \quad \text{for all } i$$

This is called **geometry optimization** or **energy minimization**. Starting from an initial guess (an experimental crystal structure, a hand-built molecule, or a random configuration), the algorithm walks downhill on the PES until the forces are negligibly small.

**Convergence criterion:** Typically we require the maximum force component on any atom to fall below a threshold, usually $f_{\max} < 0.05$ eV/Å for moderate accuracy or $< 0.01$ eV/Å for high accuracy.

#### Common Algorithms

**Steepest Descent:** Move each atom a small step in the direction of its force (downhill). Simple and robust, but slow — it zigzags through narrow valleys instead of moving along them efficiently.

**Conjugate Gradient:** Smarter choice of step direction. Instead of the pure gradient, it uses a direction that incorporates information from previous steps to avoid redundant motion. Typically 10–50× faster than steepest descent.

**BFGS (Broyden–Fletcher–Goldfarb–Shanno):** Builds an approximate second-derivative matrix (the Hessian) from successive gradient evaluations. Finds the minimum in far fewer steps than gradient-only methods, making it the default choice in ASE. It can struggle if the initial structure is very far from the minimum, in which case steepest descent or FIRE is used first.

**FIRE (Fast Inertial Relaxation Engine):** A molecular-dynamics-based minimizer that applies velocity damping. Works well for large, disordered systems.

#### In Practice

In **ASE**, geometry optimization uses the `ase.optimize` module:

```python
from ase.optimize import BFGS

opt = BFGS(atoms, logfile='optimization.log')
opt.run(fmax=0.01)   # Converge until max force < 0.01 eV/Å
```

In **LAMMPS**, the `minimize` command handles this:

```lammps
minimize 1e-6 1e-8 10000 100000
#        etol  ftol  maxiter maxeval
```

where `etol` is the fractional energy change tolerance and `ftol` is the force tolerance in energy/distance units.

:::{important} Local vs. Global Minimum
Energy minimization finds the **nearest** local minimum, not the global one. Starting from the experimental crystal structure almost always reaches the correct local minimum. Starting from a random configuration can trap the system in a high-energy metastable state.
:::

---

### 3. Energetic Quantities as Differences

This is the central concept of the lecture: **absolute energies from interatomic potentials are meaningless**. They depend on the arbitrary choice of reference state. Only differences between energies — calculated consistently with the same potential — carry physical meaning.

The art is knowing how to set up the reference states correctly.

#### 3.1 Cohesive Energy

**Definition:** The energy released per atom when isolated gas-phase atoms condense into a crystal.

$$E_{\text{coh}} = \frac{N \cdot E_{\text{atom}} - E_{\text{crystal}}}{N}$$

where $E_{\text{atom}}$ is the energy of a single isolated atom (computed by placing one atom in a large empty box) and $E_{\text{crystal}}$ is the total energy of the $N$-atom crystal after geometry optimization.

A positive $E_{\text{coh}}$ means the crystal is more stable than isolated atoms. For copper, $E_{\text{coh}} \approx 3.49$ eV/atom experimentally; for tungsten, $\approx 8.9$ eV/atom — reflecting tungsten's much stronger metallic bonding and correspondingly high melting point.

#### 3.2 Formation Energy

**Definition:** The energy change when a compound forms from its constituent elements in their standard reference structures.

For a binary compound $\text{A}_x\text{B}_y$:

$$\Delta E_f = \frac{E(\text{A}_x\text{B}_y) - x \cdot E(\text{A}_{\text{bulk}}) - y \cdot E(\text{B}_{\text{bulk}})}{x + y}$$

A **negative** $\Delta E_f$ means the compound is thermodynamically stable relative to a mechanical mixture of the pure elements. Formation energies are the backbone of computational thermodynamics — they tell you whether a new material will spontaneously form or decompose.

#### 3.3 Binding / Interaction Energy

**Definition:** The energy released when two subsystems A and B come together.

$$E_{\text{bind}} = E_{\text{AB}} - E_{\text{A}} - E_{\text{B}}$$

Negative $E_{\text{bind}}$ means binding is favorable. Applications include molecular adsorption on surfaces, dimerization of molecules, and coating adhesion.

#### 3.4 Reaction Energy

**Definition:** The energy change of a chemical reaction.

$$\Delta E_{\text{rxn}} = \sum_{\text{products}} E_i - \sum_{\text{reactants}} E_j$$

For H$_2$ dissociation:

$$\text{H}_2 \rightarrow 2\text{H}, \quad \Delta E = 2 E(\text{H}) - E(\text{H}_2) \approx +4.52 \text{ eV}$$

The positive sign means this reaction is endothermic at 0 K — energy must be supplied to break the H–H bond. As we will see in Part 2, the temperature-dependent answer is more nuanced.

---

### 4. Phase Stability at 0 K

Many materials can exist in different crystal structures (polymorphs) depending on conditions. The most stable structure at a given condition is the one with the lowest free energy. At 0 K this reduces to the lowest total energy, so we can determine phase stability with simple static energy calculations.

**Example: Iron.** Iron exists in two common structures: BCC (α-iron, stable at room temperature up to 912°C) and FCC (γ-iron, stable from 912°C to 1394°C). A static calculation predicts $E_{\text{BCC}} < E_{\text{FCC}}$, consistent with BCC being the room-temperature ground state. But this immediately raises a question: if BCC is lower in energy, why does iron transform to FCC at 912°C?

The answer is that energy alone does not determine stability at finite temperature. This is the subject of Part 2.

#### Why ML Potentials Are Well-Suited for Static Energetics

ML potentials combine near-DFT accuracy with speeds fast enough to run geometry optimizations interactively. A calculation that would take hours with DFT takes seconds with MACE. This makes it practical to compare many structures, compositions, or configurations in a single session — the kind of systematic exploration needed for alloy design and phase diagram construction.

---

## Part 2: Finite-Temperature Free Energy

### 5. Why 0 K Energies Are Incomplete

The static energy comparison above told us that BCC iron is more stable than FCC iron. Yet at 912°C iron transforms to FCC. How can a material prefer a higher-energy structure?

At finite temperature, atoms are not frozen at their energy minima — they vibrate, rotate, and translate. This thermal motion generates **entropy**, and entropy changes the landscape of stability.

The correct thermodynamic potential at constant temperature and volume is the **Helmholtz free energy**:

$$F = E - TS$$

At constant temperature and pressure (more common experimentally):

$$G = H - TS = E + PV - TS$$

A phase transition occurs when $\Delta G = 0$. A phase with slightly higher energy can be thermodynamically preferred if it has sufficiently higher entropy, provided $T\Delta S > \Delta E$. This is the fundamental competition: **energy favors the tightly packed, stiff structure; entropy favors the softer, more disordered one.**

---

### 6. Contributions to Free Energy

For molecular systems in the gas phase, the dominant free energy contributions are:

**Translational entropy:** Each molecule moving freely in a volume contributes translational entropy. Reactions that produce more independent particles gain large translational entropy — this is the driving force behind the H$_2$ dissociation example below.

**Rotational entropy:** Non-spherical molecules store energy in rotation and contribute rotational entropy. Linear molecules like H$_2$ have two rotational degrees of freedom.

**Vibrational contributions:** Molecules vibrate at discrete quantum frequencies. Each vibrational mode contributes both zero-point energy (present even at $T = 0$, largest for light atoms with stiff bonds) and thermal energy (populated as $T$ increases). The H–H stretch in H$_2$ has a characteristic temperature of ~6000 K — it is not thermally populated at room temperature but contributes a significant zero-point correction to the 0 K energy.

**Electronic:** The electronic ground state energy dominates at 0 K. Electronic excitations are generally small corrections at typical temperatures.

---

### 7. Ideal Gas Thermochemistry

For isolated molecules in the gas phase, thermodynamic quantities can be computed accurately within the **harmonic approximation**. The procedure is:

1. **Optimize geometry**: Find the minimum-energy structure where all forces vanish.
2. **Compute vibrational frequencies**: Displace each atom slightly and measure how forces change (finite-difference Hessian). For a non-linear molecule with $N$ atoms there are $3N - 6$ vibrational modes; for a linear molecule, $3N - 5$.
3. **Apply statistical mechanics**: Use the quantum harmonic oscillator partition function for each vibrational mode, plus standard translational and rotational partition functions.

This gives:

$$G(T, P) = H(T) - T S(T)$$

where $H(T)$ includes electronic energy, zero-point energy, thermal vibrational energy, and $k_B T$ (the $PV$ term for an ideal gas), and $S(T)$ includes translational, rotational, and vibrational entropy. ASE's `IdealGasThermo` class implements all of this automatically once you provide the vibrational frequencies and optimized geometry.

---

### 8. H₂ Dissociation: ΔG(T)

As a concrete example, consider:

$$\text{H}_2 \rightarrow 2\text{H}$$

**At 0 K (static energy):**

$$\Delta E = 2 E(\text{H}) - E(\text{H}_2) \approx +4.52 \text{ eV}$$

The reaction is strongly endothermic. The H–H bond must be broken, costing energy.

**Zero-point correction:**

H$_2$ has a large zero-point energy (high-frequency H–H stretch at ~4160 cm$^{-1}$). Isolated H atoms have no vibrational ZPE. Including this correction gives:

$$\Delta H_0 \approx +4.27 \text{ eV}$$

**At finite temperature:**

The entropy change $\Delta S$ is large and positive — two independent H atoms have far more translational and electronic degrees of freedom than one H$_2$ molecule. The $-T\Delta S$ term favors dissociation and grows with temperature.

$$\Delta G(T) = \Delta H(T) - T \Delta S(T)$$

At 300 K, $T\Delta S \approx 0.4$ eV — not enough to overcome $\Delta H \approx 4.3$ eV, so H$_2$ is strongly stable. The crossover $\Delta G = 0$ occurs at roughly 3000–4000 K at standard pressure, consistent with H$_2$ dissociating in high-temperature plasmas but not under normal laboratory conditions.

This illustrates a general principle: **a reaction that is endothermic but produces more particles can become spontaneous at high temperature**, driven by the entropy of the products. The crossover temperature is approximately:

$$T^* \approx \frac{\Delta H}{\Delta S}$$

---

### 9. Summary

**Part 1 (0 K):** Compare static energies after geometry optimization. Simple, computationally cheap, and directly interpretable. Valid only at absolute zero, but often gives the correct qualitative ordering of stable phases.

**Part 2 (Finite T):** Include entropy contributions using ideal gas statistical mechanics (for molecules) or more advanced methods (for solids). Necessary to explain reaction spontaneity, phase transitions, and the stability of molecules at elevated temperatures.

ML potentials make both types of calculation accessible at near-DFT accuracy and interactive speed. In Tutorial 14 you will carry out each of these calculations in Python using ASE and MACE.

---

### Further Reading

- Frenkel & Smit, *Understanding Molecular Simulation*, Chapter 5 (free energy methods)
- Sholl & Steckel, *Density Functional Theory: A Practical Introduction*, Chapter 6 (thermochemistry)
- ASE documentation: [ase.eos](https://wiki.fysik.dtu.dk/ase/ase/eos.html), [ase.thermochemistry](https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html)
