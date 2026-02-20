# Tutorial 7: Thermostats — Berendsen vs. Nosé-Hoover Chains

**Objective:** Understand *why* thermostat choice matters by implementing both algorithms from scratch and directly comparing their behavior.

:::{note} Prerequisites
This tutorial builds directly on the MD engine from Tutorial 6. Run all cells from that tutorial first so that `compute_forces`, `velocity_verlet_step`, `positions`, `velocities`, `BOX_SIZE`, `SIGMA`, `EPSILON`, and `MASS` are all available in your notebook's memory.
:::

---

## 1. Why Thermostats Matter

In Tutorial 6, we built an NVE (constant **N**umber, **V**olume, **E**nergy) simulation. Newton's laws naturally conserve total energy, so the temperature fluctuates freely — exactly as it would in a perfectly insulated box.

But real laboratory experiments are almost never done in a perfectly insulated box. A test tube sitting on a bench is in contact with the room. Heat flows in and out to keep temperature roughly constant. That is an NVT (**N**umber, **V**olume, **T**emperature) simulation — and to run one, we need a **thermostat**.

The thermostat's job sounds simple: keep temperature at a target value $T_0$. But *how* it does this has deep consequences.

> **The Central Question of this Tutorial:**
> *Does the thermostat produce a physically correct distribution of microstates, or does it just make the average temperature look right?*

We will find that the Berendsen thermostat achieves the second goal but fails the first, while Nosé-Hoover chains achieve both.

---

## 2. Measuring Temperature in a Simulation

Before we can control temperature, we need to measure it. From the Equipartition Theorem (covered in the lecture notes), the instantaneous kinetic temperature of a 2D system with $N$ particles is:

$$T_{inst} = \frac{m}{N_f \cdot k_B} \sum_{i=1}^{N} v_i^2$$

where $N_f = 2N - 2$ is the number of degrees of freedom (2 per particle in 2D, minus 2 for the fixed center-of-mass). In our reduced LJ units, $k_B = 1$ and $m = 1$, which simplifies this considerably.

**Run this cell to define the temperature calculator:**

```python
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 12

# --- In reduced LJ units: kB = 1, mass = 1 ---
# Degrees of freedom for 2D: 2*N - 2 (subtract 2 for fixed CoM)
def compute_temperature(vel, N):
    """
    Compute instantaneous kinetic temperature in reduced LJ units.
    In 2D: T = sum(v^2) / (N_dof), where N_dof = 2N - 2
    """
    N_dof = 2 * N - 2
    ke = 0.5 * MASS * np.sum(vel**2)
    return (2.0 * ke) / N_dof

print("Temperature calculator defined.")
print(f"Using N={N_PARTICLES} particles → N_dof = {2*N_PARTICLES - 2}")
```

---

## 3. The Berendsen Thermostat

### 3.1 The Idea

The Berendsen thermostat is elegant in its simplicity. Every timestep, it asks: *"Is the system too hot or too cold?"* Then it gently rescales all velocities to nudge the temperature toward the target.

The rescaling factor $\lambda$ is:

$$\lambda = \sqrt{1 + \frac{\Delta t}{\tau_T} \left( \frac{T_0}{T_{inst}} - 1 \right)}$$

**Intuition:**
* If $T_{inst} > T_0$ (too hot): $\lambda < 1$, velocities are shrunk, temperature drops.
* If $T_{inst} < T_0$ (too cold): $\lambda > 1$, velocities are stretched, temperature rises.
* The coupling time $\tau_T$ controls how aggressively the bath acts.
  * Small $\tau_T$: tight control, temperature barely fluctuates.
  * Large $\tau_T$: loose control, temperature fluctuates more naturally.

### 3.2 Implementation

**Run this cell to define the Berendsen thermostat:**

```python
def berendsen_thermostat(vel, T_current, T_target, dt, tau_T):
    """
    Berendsen velocity-rescaling thermostat.

    Parameters
    ----------
    vel      : (N, 2) array of velocities
    T_current: instantaneous temperature at this step
    T_target : desired temperature
    dt       : timestep
    tau_T    : coupling time constant (larger = weaker coupling)

    Returns
    -------
    vel_scaled : rescaled velocities
    """
    if T_current < 1e-10:
        return vel   # Avoid divide-by-zero if system is frozen

    # The rescaling factor (Eq. from lecture notes)
    lambda_sq = 1.0 + (dt / tau_T) * ((T_target / T_current) - 1.0)

    # Safety clamp: prevent extreme rescaling in the first few steps
    lambda_sq = max(0.1, min(lambda_sq, 10.0))

    return vel * np.sqrt(lambda_sq)

print("Berendsen thermostat defined.")
```

### 3.3 Running a Berendsen NVT Simulation

Now let's run the simulation and record the temperature at every step.

**Run this cell:**

```python
# --- PARAMETERS ---
T_TARGET = 1.0      # Target temperature in reduced LJ units
DT       = 0.005    # Timestep
TAU_T    = 0.5      # Coupling time (in reduced time units)
STEPS    = 10000    # Number of steps

# --- INITIALIZE from Tutorial 6 starting state ---
pos_B = positions.copy()
vel_B = velocities.copy()

# Rescale initial velocities to be close to T_TARGET
T_init = compute_temperature(vel_B, N_PARTICLES)
if T_init > 1e-10:
    vel_B *= np.sqrt(T_TARGET / T_init)

acc_B, _ = compute_forces(pos_B)

# --- DATA STORAGE ---
temps_B  = []
ke_B     = []
pe_B     = []
etot_B   = []

# --- MAIN LOOP ---
print(f"Running Berendsen NVT: {STEPS} steps, T_target={T_TARGET}, tau_T={TAU_T}")

for step in range(STEPS):
    # 1. Standard Velocity Verlet step
    pos_B, vel_B, acc_B, pe = velocity_verlet_step(pos_B, vel_B, acc_B, DT)

    # 2. Apply Berendsen thermostat AFTER the full VV step
    T_inst = compute_temperature(vel_B, N_PARTICLES)
    vel_B  = berendsen_thermostat(vel_B, T_inst, T_TARGET, DT, TAU_T)

    # 3. Record thermodynamic data
    T_new  = compute_temperature(vel_B, N_PARTICLES)
    ke     = 0.5 * MASS * np.sum(vel_B**2)
    temps_B.append(T_new)
    ke_B.append(ke)
    pe_B.append(pe)
    etot_B.append(ke + pe)

temps_B = np.array(temps_B)
ke_B    = np.array(ke_B)
pe_B    = np.array(pe_B)
etot_B  = np.array(etot_B)

print(f"Done. Mean T = {np.mean(temps_B):.4f} (target: {T_TARGET})")
print(f"Std  T = {np.std(temps_B):.4f}")
```

---

## 4. The Nosé-Hoover Chain Thermostat

### 4.1 The Problem with Berendsen

The Berendsen thermostat works by directly manipulating velocities every step. This is like a hand pressing on a spring — it forces the outcome but leaves no room for the natural fluctuations that the Canonical Ensemble (NVT) requires.

In the true Canonical Ensemble, temperature **should** fluctuate. The variance of these fluctuations is not random noise — it is physically meaningful. It is related to the heat capacity of the system:

$$\langle (\Delta T)^2 \rangle = \frac{2 T_0^2}{N_f}$$

Berendsen suppresses this variance. Nosé-Hoover respects it.

### 4.2 The Core Idea

The Nosé-Hoover approach introduces a fictitious "heat bath" degree of freedom $\xi$ (the thermostat variable) that couples to the atoms. Instead of rescaling velocities externally, it adds a friction-like term to the equations of motion:

$$\dot{\mathbf{v}}_i = \frac{\mathbf{F}_i}{m} - \xi \mathbf{v}_i$$

The thermostat variable $\xi$ itself evolves according to its own equation:

$$\dot{\xi} = \frac{1}{Q} \left( \sum_i m v_i^2 - N_f k_B T_0 \right)$$

**Intuition:** The term in brackets is the difference between the current kinetic energy and its target value. If the system is too hot, $\dot{\xi}$ is positive, so $\xi$ grows, which increases the friction on the atoms, slowing them down. If the system is too cold, the reverse happens. The "mass" $Q$ controls the inertia of the bath — how quickly it responds.

**Chains:** A single Nosé-Hoover thermostat can sometimes fail to ergodically sample all states (especially for stiff or small systems). The fix is to chain multiple thermostats together, where the second thermostat controls the first, and so on. We will implement a 2-chain version, which is robust for our 2D system.

### 4.3 Implementation

**Run this cell to define the Nosé-Hoover Chain integrator:**

```python
def nose_hoover_chain_step(pos, vel, acc, xi1, xi2, p_xi1, p_xi2, Q1, Q2,
                            dt, T_target, N_dof, force_fn):
    """
    One timestep of Velocity Verlet + Nosé-Hoover Chain (2 links).

    The chain is integrated using the Liouville operator splitting (RESPA-style),
    which keeps the extended system symplectic and time-reversible.

    Parameters
    ----------
    pos, vel, acc   : positions, velocities, accelerations  (N, 2)
    xi1, xi2        : thermostat position variables (scalars)
    p_xi1, p_xi2    : thermostat momentum variables (scalars)
    Q1, Q2          : thermostat masses (scalars)
    dt              : timestep
    T_target        : target temperature (kB = 1 in reduced units)
    N_dof           : degrees of freedom (= 2N - 2 in 2D)
    force_fn        : function(pos) -> (acc, pe)

    Returns
    -------
    Updated pos, vel, acc, xi1, xi2, p_xi1, p_xi2, pe
    """
    # --- 1. Half-step update of chain momenta (outermost first) ---
    # Update p_xi2 using friction from p_xi1
    KE2 = p_xi1**2 / (2.0 * Q1)
    p_xi2 += 0.5 * dt * (KE2 - T_target)

    # Update p_xi1 using friction from atoms and p_xi2
    KE_atoms = 0.5 * MASS * np.sum(vel**2)
    scale1 = np.exp(-0.5 * dt * p_xi2 / Q2)
    p_xi1  = p_xi1 * scale1 + 0.5 * dt * scale1 * (2.0 * KE_atoms - N_dof * T_target)
    p_xi1 *= scale1     # symmetric: apply scale factor twice (Yoshida splitting)

    # --- 2. Rescale velocities using the chain friction ---
    vel_scale = np.exp(-0.5 * dt * p_xi1 / Q1)
    vel = vel * vel_scale

    # --- 3. Standard Velocity Verlet: half-kick → drift → force → half-kick ---
    vel_half = vel + 0.5 * acc * dt
    pos_new  = (pos + vel_half * dt) % BOX_SIZE
    acc_new, pe = force_fn(pos_new)
    vel_new  = vel_half + 0.5 * acc_new * dt

    # --- 4. Half-step update of chain momenta (reverse order) ---
    # Rescale again (second half)
    vel_new = vel_new * vel_scale

    KE_atoms_new = 0.5 * MASS * np.sum(vel_new**2)
    scale1_rev   = np.exp(-0.5 * dt * p_xi2 / Q2)
    p_xi1 = p_xi1 * scale1_rev + 0.5 * dt * scale1_rev * (2.0 * KE_atoms_new - N_dof * T_target)
    p_xi1 *= scale1_rev

    KE2_new = p_xi1**2 / (2.0 * Q1)
    p_xi2  += 0.5 * dt * (KE2_new - T_target)

    # --- 5. Advance thermostat positions ---
    xi1 += dt * p_xi1 / Q1
    xi2 += dt * p_xi2 / Q2

    return pos_new, vel_new, acc_new, xi1, xi2, p_xi1, p_xi2, pe

print("Nosé-Hoover Chain integrator defined.")
```

### 4.4 Running a Nosé-Hoover NVT Simulation

**Run this cell:**

```python
# --- PARAMETERS (same physics, same target T) ---
T_TARGET = 1.0
DT       = 0.005
STEPS    = 10000

# Thermostat "mass" Q: controls the inertia of the heat bath.
# Rule of thumb: Q ≈ N_dof * kB * T_target * tau_T^2
# We use tau_T = 0.2 (tighter coupling for better ergodic sampling)
N_DOF   = 2 * N_PARTICLES - 2
TAU_NH  = 0.2
Q1 = N_DOF * T_TARGET * TAU_NH**2   # Mass of chain link 1
Q2 = T_TARGET * TAU_NH**2           # Mass of chain link 2

print(f"NH Chain masses: Q1={Q1:.2f}, Q2={Q2:.4f}")

# --- INITIALIZE (same starting state for a fair comparison) ---
pos_NH  = positions.copy()
vel_NH  = velocities.copy()

# Rescale to T_TARGET
T_init = compute_temperature(vel_NH, N_PARTICLES)
if T_init > 1e-10:
    vel_NH *= np.sqrt(T_TARGET / T_init)

acc_NH, _ = compute_forces(pos_NH)

# Initialize chain variables at rest
xi1, xi2     = 0.0, 0.0
p_xi1, p_xi2 = 0.0, 0.0

# --- DATA STORAGE ---
temps_NH = []
ke_NH    = []
pe_NH    = []
etot_NH  = []

# --- MAIN LOOP ---
print(f"Running Nosé-Hoover Chain NVT: {STEPS} steps, T_target={T_TARGET}")

for step in range(STEPS):
    pos_NH, vel_NH, acc_NH, xi1, xi2, p_xi1, p_xi2, pe = nose_hoover_chain_step(
        pos_NH, vel_NH, acc_NH, xi1, xi2, p_xi1, p_xi2,
        Q1, Q2, DT, T_TARGET, N_DOF, compute_forces
    )

    T_inst = compute_temperature(vel_NH, N_PARTICLES)
    ke     = 0.5 * MASS * np.sum(vel_NH**2)
    temps_NH.append(T_inst)
    ke_NH.append(ke)
    pe_NH.append(pe)
    etot_NH.append(ke + pe)

temps_NH = np.array(temps_NH)
ke_NH    = np.array(ke_NH)
pe_NH    = np.array(pe_NH)
etot_NH  = np.array(etot_NH)

print(f"Done. Mean T = {np.mean(temps_NH):.4f} (target: {T_TARGET})")
print(f"Std  T = {np.std(temps_NH):.4f}")
```

---

## 5. Head-to-Head Comparison

Now we compare the two thermostats side by side. Run the cells in order.

### 5.1 Temperature Traces

The most immediate question: does each thermostat maintain the target temperature?

**Run this cell:**

```python
steps_B  = np.arange(len(temps_B))
steps_NH = np.arange(len(temps_NH))

fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=False)

# --- Berendsen ---
axes[0].plot(steps_B, temps_B, color='steelblue', alpha=0.7, lw=0.8, label='Instantaneous T')
axes[0].axhline(T_TARGET,         color='red',  lw=1.5, ls='--', label=f'Target ({T_TARGET})')
axes[0].axhline(np.mean(temps_B), color='navy', lw=1.5, ls='-',
                label=f'Mean = {np.mean(temps_B):.3f}')
axes[0].fill_between(steps_B,
                     np.mean(temps_B) - np.std(temps_B),
                     np.mean(temps_B) + np.std(temps_B),
                     alpha=0.15, color='steelblue', label=f'±1σ (σ={np.std(temps_B):.4f})')
axes[0].set_ylabel('Temperature (LJ units)', fontsize=11)
axes[0].set_title('Berendsen Thermostat — Temperature Trace', fontsize=13)
axes[0].legend(loc='upper right', fontsize=9)
axes[0].set_ylim(0, 2.5)

# --- Nosé-Hoover ---
axes[1].plot(steps_NH, temps_NH, color='darkorange', alpha=0.7, lw=0.8, label='Instantaneous T')
axes[1].axhline(T_TARGET,          color='red',         lw=1.5, ls='--', label=f'Target ({T_TARGET})')
axes[1].axhline(np.mean(temps_NH), color='saddlebrown', lw=1.5, ls='-',
                label=f'Mean = {np.mean(temps_NH):.3f}')
axes[1].fill_between(steps_NH,
                     np.mean(temps_NH) - np.std(temps_NH),
                     np.mean(temps_NH) + np.std(temps_NH),
                     alpha=0.15, color='darkorange', label=f'±1σ (σ={np.std(temps_NH):.4f})')
axes[1].set_ylabel('Temperature (LJ units)', fontsize=11)
axes[1].set_xlabel('Step', fontsize=11)
axes[1].set_title('Nosé-Hoover Chain — Temperature Trace', fontsize=13)
axes[1].legend(loc='upper right', fontsize=9)
axes[1].set_ylim(0, 2.5)

plt.suptitle('Thermostat Comparison: Temperature Control', fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
plt.show()

print("\nObservation:")
print(f"  Berendsen   std(T) = {np.std(temps_B):.4f}")
print(f"  Nosé-Hoover std(T) = {np.std(temps_NH):.4f}")
print(f"\n  Berendsen suppresses fluctuations — its σ(T) is artificially small.")
print(f"  Nosé-Hoover allows natural fluctuations expected from statistical mechanics.")
```

### 5.2 The Canonical Ensemble Test: Kinetic Energy Distribution

This is the most important diagnostic. The Canonical Ensemble (NVT) predicts that the **kinetic energy** of a system with $N_f$ degrees of freedom follows a **chi-squared distribution** with $N_f$ degrees of freedom:

$$P(K) \propto K^{N_f/2 - 1} \, e^{-K / k_B T_0}$$

This is a specific, testable prediction. Berendsen does not satisfy it. Nosé-Hoover chains do.

**Run this cell:**

```python
from scipy.stats import chi2

# --- Theoretical chi-squared distribution ---
# For N_dof degrees of freedom, KE ~ (kB*T/2) * chi2(N_dof)
# So KE / (kB*T/2) ~ chi2(N_dof)
N_dof  = 2 * N_PARTICLES - 2
scale  = T_TARGET / 2.0          # kB*T/2 in reduced units

# Discard first 1000 steps (equilibration)
ke_B_eq  = ke_B[1000:]
ke_NH_eq = ke_NH[1000:]

ke_range = np.linspace(0, max(ke_B_eq.max(), ke_NH_eq.max()) * 1.1, 300)
theory   = chi2.pdf(ke_range / scale, df=N_dof) / scale   # Proper normalization

fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

# --- Berendsen ---
axes[0].hist(ke_B_eq, bins=50, density=True, color='steelblue', alpha=0.6,
             edgecolor='black', lw=0.3, label='Berendsen')
axes[0].plot(ke_range, theory, 'r-', lw=2.5, label=f'Canonical (χ² df={N_dof})')
axes[0].set_xlabel('Kinetic Energy (LJ units)', fontsize=11)
axes[0].set_ylabel('Probability Density', fontsize=11)
axes[0].set_title('Berendsen: KE Distribution', fontsize=13)
axes[0].legend(fontsize=10)

# --- Nosé-Hoover ---
axes[1].hist(ke_NH_eq, bins=50, density=True, color='darkorange', alpha=0.6,
             edgecolor='black', lw=0.3, label='Nosé-Hoover Chain')
axes[1].plot(ke_range, theory, 'r-', lw=2.5, label=f'Canonical (χ² df={N_dof})')
axes[1].set_xlabel('Kinetic Energy (LJ units)', fontsize=11)
axes[1].set_title('Nosé-Hoover Chain: KE Distribution', fontsize=13)
axes[1].legend(fontsize=10)

plt.suptitle('Canonical Ensemble Test: Does the KE follow χ²?', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.show()

print("\nInterpretation:")
print("  If the histogram matches the red curve → correct canonical sampling.")
print("  Berendsen: narrow histogram, does NOT match the red curve.")
print("  Nosé-Hoover: broader histogram, matches the red curve.")
```

### 5.3 Energy Conservation Diagnostic

In an NVE simulation, total energy ($K + U$) is conserved. With a thermostat, the system exchanges energy with the bath, so total energy should *fluctuate*, not be conserved. However, the Nosé-Hoover chain has a **conserved extended energy** that includes the bath. We can check this as a numerical sanity test.

The Nosé-Hoover conserved quantity is:
$$H_{ext} = K + U + \frac{p_{\xi_1}^2}{2Q_1} + \frac{p_{\xi_2}^2}{2Q_2} + N_f k_B T_0 \xi_1 + k_B T_0 \xi_2$$

For Berendsen, there is no such conserved quantity — energy is simply injected or removed at every step. This is why it cannot produce the correct ensemble.

**Run this cell:**

```python
def running_average(arr, window=200):
    """Simple moving average using convolution."""
    kernel = np.ones(window) / window
    return np.convolve(arr, kernel, mode='valid')

steps_B  = np.arange(len(etot_B))
steps_NH = np.arange(len(etot_NH))

window = 200

ra_B  = running_average(etot_B,  window)
ra_NH = running_average(etot_NH, window)

# x-axes for the running averages (shorter due to convolution 'valid' mode)
steps_ra_B  = steps_B[window - 1:]
steps_ra_NH = steps_NH[window - 1:]

fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=False)

# --- Berendsen ---
axes[0].plot(steps_B, etot_B, color='steelblue', lw=0.6, alpha=0.4, label='Instantaneous K+U')
axes[0].plot(steps_ra_B, ra_B, color='navy', lw=2.0, label=f'Running avg (window={window})')
axes[0].set_ylabel('K + U (LJ units)', fontsize=11)
axes[0].set_title('Berendsen: Total Mechanical Energy — running average drifts freely', fontsize=11)
axes[0].legend(fontsize=9)

# --- Nosé-Hoover ---
axes[1].plot(steps_NH, etot_NH, color='darkorange', lw=0.6, alpha=0.4, label='Instantaneous K+U')
axes[1].plot(steps_ra_NH, ra_NH, color='saddlebrown', lw=2.0, label=f'Running avg (window={window})')
axes[1].set_ylabel('K + U (LJ units)', fontsize=11)
axes[1].set_xlabel('Step', fontsize=11)
axes[1].set_title('Nosé-Hoover Chain: Total Mechanical Energy — fluctuates around a stable mean', fontsize=11)
axes[1].legend(fontsize=9)

plt.suptitle('Mechanical Energy: Berendsen vs. Nosé-Hoover', fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
plt.show()

print("\nObservation:")
print(f"  Berendsen   — mean K+U = {np.mean(etot_B):.3f},  std = {np.std(etot_B):.4f}")
print(f"  Nosé-Hoover — mean K+U = {np.mean(etot_NH):.3f},  std = {np.std(etot_NH):.4f}")
print(f"\n  Berendsen has no conserved quantity — the bath injects/removes energy ad hoc.")
print(f"  Nosé-Hoover fluctuates but has a conserved extended Hamiltonian (bath included).")
```

---

## 6. Exercise: Effect of Coupling Strength

The coupling time $\tau_T$ (Berendsen) or equivalently the choice of $Q$ (Nosé-Hoover) is the most important parameter to tune.

**Task:** Run the cell below. It runs Berendsen at three different $\tau_T$ values and plots the temperature traces. Observe how the thermostat behavior changes.

**Questions to answer:**
1. Which $\tau_T$ gives the tightest control? Which gives the most natural-looking fluctuations?
2. For a $\tau_T$ close to $\Delta t$, the thermostat is nearly identical to direct velocity rescaling. What happens?
3. For a very large $\tau_T$, does the mean temperature still converge to the target within 10,000 steps?

```python
# --- EXERCISE: Scan tau_T values for Berendsen ---
tau_values  = [0.05, 0.5, 5.0]
colors      = ['purple', 'steelblue', 'teal']
DT_EX       = 0.005
STEPS_EX    = 10000
T_TARGET_EX = 1.0

fig, axes = plt.subplots(len(tau_values), 1, figsize=(12, 9), sharex=True, sharey=True)

for idx, (tau, color) in enumerate(zip(tau_values, colors)):
    # Fresh start each time
    pos_ex = positions.copy()
    vel_ex = velocities.copy()
    T_i = compute_temperature(vel_ex, N_PARTICLES)
    if T_i > 1e-10:
        vel_ex *= np.sqrt(T_TARGET_EX / T_i)
    acc_ex, _ = compute_forces(pos_ex)
    temps_ex  = []

    for step in range(STEPS_EX):
        pos_ex, vel_ex, acc_ex, _ = velocity_verlet_step(pos_ex, vel_ex, acc_ex, DT_EX)
        T_inst  = compute_temperature(vel_ex, N_PARTICLES)
        vel_ex  = berendsen_thermostat(vel_ex, T_inst, T_TARGET_EX, DT_EX, tau)
        temps_ex.append(compute_temperature(vel_ex, N_PARTICLES))

    temps_ex = np.array(temps_ex)

    axes[idx].plot(temps_ex, color=color, lw=0.8, alpha=0.8)
    axes[idx].axhline(T_TARGET_EX, color='red', lw=1.5, ls='--', label='Target')
    axes[idx].set_ylabel('T (LJ units)', fontsize=10)
    axes[idx].set_title(
        f'τ_T = {tau}  →  mean T = {np.mean(temps_ex):.3f}, σ(T) = {np.std(temps_ex):.4f}',
        fontsize=11
    )
    axes[idx].set_ylim(-0.1, 3.0)
    axes[idx].legend(fontsize=9)

axes[-1].set_xlabel('Step', fontsize=11)
plt.suptitle('Berendsen Thermostat: Effect of Coupling Time τ_T', fontsize=13, fontweight='bold')
plt.tight_layout()
plt.show()
```

---

## 7. Summary: When to Use Each Thermostat

The table below summarizes the practical guidance from what we observed in this tutorial. This directly mirrors the guidance in the lecture notes.

| Property | Berendsen | Nosé-Hoover Chain |
| :--- | :--- | :--- |
| **Mean temperature** | ✅ Converges to $T_0$ | ✅ Converges to $T_0$ |
| **Temperature fluctuations** | ❌ Artificially suppressed | ✅ Physically correct |
| **KE distribution** | ❌ Too narrow (not χ²) | ✅ Correct canonical (χ²) |
| **Conserved quantity** | ❌ None (energy injected ad hoc) | ✅ Extended Hamiltonian conserved |
| **Free energy calculations** | ❌ **Do not use** | ✅ Correct |
| **Equilibration** | ✅ Fast, robust, hard to break | ⚠️ Slower equilibration |
| **Small/stiff systems** | ✅ Reliable | ⚠️ May need more chain links |
| **Typical use** | Heating/equilibration phase | Production data collection |

:::{important} The Practical Workflow
The typical protocol (used in LAMMPS, GROMACS, NAMD) is:
1. **Energy Minimization** — remove clashes.
2. **NVT with Berendsen** — heat the system quickly to the target temperature.
3. **NPT with Berendsen** — equilibrate the density.
4. **Switch to Nosé-Hoover** — run the production simulation from which you extract scientific results.

Never use Berendsen data for thermodynamic properties like free energy, heat capacity, or entropy.
:::

:::{note} Further Reading
* Berendsen et al., *J. Chem. Phys.* **81**, 3684 (1984) — the original paper.
* Nosé, *Mol. Phys.* **52**, 255 (1984) and Hoover, *Phys. Rev. A* **31**, 1695 (1985) — the foundational papers.
* Martyna, Klein & Tuckerman, *J. Chem. Phys.* **97**, 2635 (1992) — the chains extension.
* Frenkel & Smit, *Understanding Molecular Simulation*, Chapter 6 — the clearest textbook treatment.
:::
