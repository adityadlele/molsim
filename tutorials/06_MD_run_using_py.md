
# Tutorial 6: Building an MD Engine from Scratch

**Objective:** Understand the algorithms inside LAMMPS by building a modular MD engine in Python, step-by-step.

## 1. Introduction
In previous tutorials, we treated LAMMPS as a "black box." Now, we will open that box.

We will write a Python script that explicitly performs the three core tasks of any MD simulation:
1.  **Initialize:** Create atoms with positions and velocities.
2.  **Calculate Forces:** Use the Lennard-Jones potential.
3.  **Integrate:** Move atoms using Newton's laws (Velocity Verlet).

---

## Step 1: Imports & Global Settings
First, we need to load the necessary libraries and define the "Universe" constants.
* **Box Size:** The width of our periodic container.
* **Sigma/Epsilon:** The fundamental units of size and energy for our atoms.

**Run this cell first:**
```python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython.display import HTML

# Make the plots look professional
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.figsize'] = (6, 6)
plt.rcParams['font.size'] = 12

# GLOBAL CONSTANTS (The "Universe" Settings)
BOX_SIZE = 10.0      # The size of our container
SIGMA = 1.0          # The size of an atom
EPSILON = 1.0        # The strength of attraction
MASS = 1.0           # The mass of an atom

```

---

## Step 2: The Initial State

An MD simulation needs to know two things about every atom at :

1. **Position ():** Where is it?
2. **Velocity ():** Where is it going?

We will generate random positions and velocities. Crucially, we must subtract the "Center of Mass" velocity so our entire gas cloud doesn't drift sideways.

**Run this cell to visualize the starting state:**

```python
# --- LOGIC: STRICT INITIALIZATION ---
N_PARTICLES = 20
MIN_DIST = 1.12 * SIGMA  # Distance where atoms repel (2^(1/6))
MAX_ATTEMPTS = 5000      

positions = np.zeros((N_PARTICLES, 2))

# 1. Place Atoms (Checking for Overlaps + Boundaries)
for i in range(N_PARTICLES):
    placed = False
    attempts = 0
    while not placed and attempts < MAX_ATTEMPTS:
        candidate = np.random.rand(2) * BOX_SIZE
        
        if i == 0:
            positions[i] = candidate
            placed = True
        else:
            # Calculate distance to all existing atoms
            delta = positions[:i] - candidate
            
            # --- PERIODIC BOUNDARY CHECK ---
            # If atoms are on opposite edges (0 and 10), they are actually close!
            delta -= np.round(delta / BOX_SIZE) * BOX_SIZE 
            
            dist_sq = np.sum(delta**2, axis=1)
            
            # If ALL existing atoms are far enough away, we accept this spot
            if np.all(dist_sq > MIN_DIST**2):
                positions[i] = candidate
                placed = True     
        attempts += 1
    
    if not placed:
        print(f"Warning: Failed to place particle {i}. Box is too full.")

# 2. Generate Velocities
# Random values between -0.5 and 0.5
velocities = (np.random.rand(N_PARTICLES, 2) - 0.5)

# Center of Mass Correction
# We subtract the average so the gas doesn't drift
velocities -= np.mean(velocities, axis=0)

# --- VISUALIZATION ---
fig, ax = plt.subplots(figsize=(6,6))
ax.set_xlim(0, BOX_SIZE); ax.set_ylim(0, BOX_SIZE)
ax.set_title(f"Initial State (N={N_PARTICLES})")
ax.set_xlabel("Position X"); ax.set_ylabel("Position Y")

# 1. Draw Atoms (Blue Circles)
ax.scatter(positions[:,0], positions[:,1], s=300, color='skyblue', edgecolors='black', alpha=0.8, label='Atom')

# 2. Draw Velocities (Red Arrows)
# Quiver plots arrows at (x,y) with components (u,v)
ax.quiver(positions[:,0], positions[:,1], 
          velocities[:,0], velocities[:,1], 
          color='red', width=0.005, scale=5, label='Velocity Vector')

ax.legend(loc='upper right')
plt.grid(True, linestyle='--', alpha=0.6)
plt.show()

```

---

## Step 3: The Forces (Lennard-Jones)

Atoms interact via forces. We use the **Lennard-Jones Potential**.

* **Too Close?** Strong Repulsion (Positive force).
* **Just Right?** Weak Attraction (Negative force).
* **Too Far?** Zero force (Cutoff).

This block defines the `compute_forces` function. It calculates the force on every pair of atoms. It also implements **Periodic Boundary Conditions** (Minimum Image Convention), ensuring atoms interact across the box edges.

**Run this cell to define the physics and see the force curve:**

```python
# --- LOGIC ---
# Pre-calculate the potential at the cutoff distance to shift the curve
CUTOFF = 3.0 * SIGMA
inv_rc2 = (SIGMA**2) / (CUTOFF**2)
inv_rc6 = inv_rc2**3
E_CUTOFF = 4 * EPSILON * (inv_rc6**2 - inv_rc6)

def compute_forces(pos):
    """
    Input: Positions of all atoms
    Output: Acceleration (Force/Mass) of all atoms, Potential Energy
    """
    N = pos.shape[0]
    forces = np.zeros_like(pos)
    potential_energy = 0.0
    
    for i in range(N):
        for j in range(i + 1, N): # Loop over pairs
            delta = pos[i] - pos[j]
            delta -= np.round(delta / BOX_SIZE) * BOX_SIZE # Periodic BC
            r_sq = np.sum(delta**2)
            
            if r_sq < 0.01: r_sq = 0.01
            
            # Optimization Cutoff
            if r_sq < CUTOFF**2:
                inv_r2 = (SIGMA**2) / r_sq
                inv_r6 = inv_r2**3
                
                # Force (unchanged by the energy shift)
                f_scalar = (24 * EPSILON / r_sq) * (2 * inv_r6**2 - inv_r6)
                
                force_vec = f_scalar * delta
                forces[i] += force_vec
                forces[j] -= force_vec
                
                # ENERGY CORRECTION: Subtract E_CUTOFF
                # This ensures Potential Energy is exactly 0.0 at the cutoff line
                pair_energy = 4 * EPSILON * (inv_r6**2 - inv_r6) - E_CUTOFF
                potential_energy += pair_energy
                
    return forces / MASS, potential_energy

# --- VISUALIZATION OF THE RULE ---
# We plot the "Shifted" potential to show it hits exactly zero
r = np.linspace(0.85, 3.2, 100)
E = 4 * EPSILON * ((SIGMA/r)**12 - (SIGMA/r)**6)
# Apply shift to plotting data too for consistency
E[r > CUTOFF] = 0
E[r <= CUTOFF] -= E_CUTOFF

plt.figure(figsize=(8,4))
plt.plot(r, E, linewidth=3, color='purple')
plt.axhline(0, color='k', linestyle='--')
plt.axvline(1.0, color='green', linestyle='--', label='Sigma')
plt.axvline(CUTOFF, color='red', linestyle=':', label='Cutoff (Shifted to 0)')
plt.title("Block 2: The Rules (Shifted Lennard-Jones)")
plt.xlabel("Distance"); plt.ylabel("Energy")
plt.legend()
plt.show()

```

---

## Step 4: The Integrator (Velocity Verlet)

We need an engine to move time forward. We use the **Velocity Verlet** algorithm.
It is more accurate than simple Euler integration () because it calculates velocity at "half-steps," which helps conserve energy.

**Run this cell to load the integrator:**

```python
# --- LOGIC ---
def velocity_verlet_step(pos, vel, acc, dt):
    """
    Moves the simulation forward by one step 'dt'
    """
    # 1. First Half-Kick: Update velocity halfway
    vel_half = vel + 0.5 * acc * dt
    
    # 2. Drift: Update position using half-velocity
    pos_new = pos + vel_half * dt
    pos_new = pos_new % BOX_SIZE # Wrap around box
    
    # 3. Calculate new forces at new positions
    acc_new, pe = compute_forces(pos_new)
    
    # 4. Second Half-Kick: Finish updating velocity
    vel_new = vel_half + 0.5 * acc_new * dt
    
    return pos_new, vel_new, acc_new, pe

print("Block 3: Integrator Loaded. The Engine is ready.")
```

---

## Step 5: The Simulation Loop

Now we combine the blocks into a loop.

1. Initialize positions.
2. Calculate initial forces.
3. Loop 200 times, calling `velocity_verlet_step` each time.
4. Plot the positions and the Energy.

**Run this cell to watch the movie:**

```python
# --- SETUP ---
DT = 0.005             # Reduced slightly for better stability
STEPS = 200            # How long to run

# CRITICAL FIX: Use the safe variables from Cell 2!
# We use .copy() so if you run this cell twice, it doesn't mess up the original data.
pos = positions.copy()
vel = velocities.copy()

# Calculate initial forces based on these SAFE positions
acc, pot_e = compute_forces(pos) 

# --- ANIMATION SETUP ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# 1. Particle View
ax1.set_xlim(0, BOX_SIZE); ax1.set_ylim(0, BOX_SIZE)
ax1.set_title(f"Simulation (N={len(pos)})")
particles = ax1.scatter([], [], s=150, c=[], cmap='coolwarm', edgecolors='k')

# 2. Energy View
ax2.set_xlim(0, STEPS)
ax2.set_title("Energy Conservation")
ax2.set_xlabel("Time Step"); ax2.set_ylabel("Energy")

# Initialize lines
line_ke, = ax2.plot([], [], 'r-', label='Kinetic')
line_pe, = ax2.plot([], [], 'b-', label='Potential')
line_tot, = ax2.plot([], [], 'k--', label='Total', lw=2)
ax2.legend()

history_ke, history_pe, history_tot = [], [], []

# Pre-calculate initial energy to start the graph at t=0
initial_ke = 0.5 * MASS * np.sum(np.sum(vel**2, axis=1))
history_ke.append(initial_ke)
history_pe.append(pot_e)
history_tot.append(initial_ke + pot_e)

def update(frame):
    global pos, vel, acc
    
    # CALL THE INTEGRATOR
    pos, vel, acc, pe = velocity_verlet_step(pos, vel, acc, DT)
    
    # CALCULATE KINETIC ENERGY
    ke = 0.5 * MASS * np.sum(np.sum(vel**2, axis=1))
    
    # STORE DATA
    history_ke.append(ke); history_pe.append(pe); history_tot.append(ke+pe)
    
    # UPDATE VISUALS
    particles.set_offsets(pos)
    particles.set_array(np.linalg.norm(vel, axis=1)) # Color by speed
    
    x = range(len(history_ke))
    line_ke.set_data(x, history_ke)
    line_pe.set_data(x, history_pe)
    line_tot.set_data(x, history_tot)
    
    # Auto-scale the energy plot
    if len(history_ke) > 1:
        min_y = min(min(history_pe), min(history_ke))
        max_y = max(max(history_pe), max(history_ke))
        ax2.set_ylim(min_y - 5, max_y + 5)
        
    return particles, line_ke, line_pe, line_tot

anim = animation.FuncAnimation(fig, update, frames=STEPS, interval=30, blit=False)
plt.close()
HTML(anim.to_jshtml())

```

---

## Exercise 1: The Broken Physics (Explosion)

What happens if we take a "Time Step" that is too large?
In this exercise, we set `DT = 0.05`.

**Task:** Run the cell below. Look at the Black Line (Total Energy).
**Observation:** Energy is not conserved! It explodes upwards. This demonstrates why MD requires small timesteps (usually 1 femtosecond).

```python
# --- EXERCISE 1 ---
DT_BROKEN = 0.05       # <--- TOO LARGE! 
STEPS = 100

# Reset
pos = np.random.rand(30, 2) * BOX_SIZE
vel = (np.random.rand(30, 2) - 0.5) * 2.0
acc, _ = compute_forces(pos)

fig, ax = plt.subplots(figsize=(6,4))
ax.set_title("Exercise 1: Broken Physics (Energy Explosion)")
ax.set_xlim(0, STEPS)
l_tot, = ax.plot([], [], 'k--', label='Total Energy', lw=2)
ax.legend()
hist_tot = []

def update_ex1(frame):
    global pos, vel, acc
    pos, vel, acc, pe = velocity_verlet_step(pos, vel, acc, DT_BROKEN)
    ke = 0.5 * MASS * np.sum(vel**2)
    hist_tot.append(ke + pe)
    l_tot.set_data(range(len(hist_tot)), hist_tot)
    ax.set_ylim(min(hist_tot)-10, max(hist_tot)+10)
    return l_tot,

anim = animation.FuncAnimation(fig, update_ex1, frames=STEPS, interval=30, blit=True)
plt.close()
HTML(anim.to_jshtml())


```

---

## Exercise 2: Solid Phase (Crystallization)

Can simple spheres form a solid crystal?
Here, we increase the density (`N=64` in a small box) and lower the temperature (very small velocities).

**Task:** Run the cell and watch the atoms. They should jiggle in place, forming a lattice.

```python
# --- EXERCISE 2 ---
N_SOLID = 64
BOX_SOLID = 8.0
DT = 0.01

# Initialize in a grid to help crystallization
grid_n = int(np.ceil(np.sqrt(N_SOLID)))
spacing = BOX_SOLID / grid_n
x = np.linspace(spacing/2, BOX_SOLID-spacing/2, grid_n)
y = np.linspace(spacing/2, BOX_SOLID-spacing/2, grid_n)
xv, yv = np.meshgrid(x, y)
pos = np.column_stack((xv.ravel(), yv.ravel()))[:N_SOLID]

# Very small velocities (Cold)
vel = (np.random.rand(N_SOLID, 2) - 0.5) * 0.1 
acc, _ = compute_forces(pos)

fig, ax = plt.subplots(figsize=(6,6))
ax.set_xlim(0, BOX_SOLID); ax.set_ylim(0, BOX_SOLID)
ax.set_title("Exercise 2: Solid Phase Formation")
parts = ax.scatter([], [], s=200, c='blue', edgecolors='k')

def update_solid(frame):
    global pos, vel, acc
    pos, vel, acc, _ = velocity_verlet_step(pos, vel, acc, DT)
    parts.set_offsets(pos)
    return parts,

anim = animation.FuncAnimation(fig, update_solid, frames=200, interval=50, blit=True)
plt.close()
HTML(anim.to_jshtml())
```

---

## Exercise 3: How Long Does MD Actually Take?

In the lecture notes, we emphasized that real simulations require **billions of timesteps** and are impossible to run on a laptop. But *why exactly*? Let's measure it ourselves.

This exercise does two things:
1. **Times** the main simulation loop for our 20-particle system.
2. **Scales** the particle count and plots how runtime grows — revealing the $O(N^2)$ bottleneck in `compute_forces`.

### Part A: Timing the Main Loop

**Run this cell to time 200 steps of the simulation:**

```python
import time

# --- SETUP (same as Step 5) ---
DT = 0.005
STEPS = 200

pos = positions.copy()
vel = velocities.copy()
acc, _ = compute_forces(pos)

# --- TIMED LOOP ---
start_time = time.perf_counter()

for step in range(STEPS):
    pos, vel, acc, _ = velocity_verlet_step(pos, vel, acc, DT)

end_time = time.perf_counter()

wall_time_s = end_time - start_time
time_per_step_ms = (wall_time_s / STEPS) * 1000

print("=" * 45)
print(f"  Particles (N)      : {N_PARTICLES}")
print(f"  Steps completed    : {STEPS}")
print(f"  Total wall time    : {wall_time_s:.3f} seconds")
print(f"  Time per step      : {time_per_step_ms:.3f} ms")
print("=" * 45)

# --- SCALE UP: How long for a REAL simulation? ---
# A real MD run targets ~1 microsecond = 1,000,000 steps at dt=1fs
REAL_STEPS = 1_000_000
projected_s = (wall_time_s / STEPS) * REAL_STEPS
projected_days = projected_s / 86400

print(f"\n  Projected time for {REAL_STEPS:,} steps:")
print(f"  → {projected_s:,.1f} seconds  ({projected_days:.1f} days)")
print(f"\n  This is why we need supercomputers!")
```

### Part B: The $O(N^2)$ Scaling Problem

The `compute_forces` function loops over every **pair** of atoms. For $N$ atoms, there are $N(N-1)/2$ pairs. Double the atoms, and you quadruple the work.

Let's measure this directly:

```python
import time
import numpy as np
import matplotlib.pyplot as plt

# --- PARAMETERS ---
# Test a range of particle counts
particle_counts = [5, 10, 20, 40, 80, 160]
TIMING_STEPS = 50          # Steps per timing run (keep short for large N)
timings = []

print(f"{'N':>6} | {'Time/step (ms)':>16} | {'Pairs':>10}")
print("-" * 38)

for N in particle_counts:
    # Place N atoms randomly (no overlap check for speed)
    test_pos = np.random.rand(N, 2) * BOX_SIZE
    test_vel = (np.random.rand(N, 2) - 0.5)
    test_acc, _ = compute_forces(test_pos)
    
    # Time TIMING_STEPS steps
    t0 = time.perf_counter()
    for _ in range(TIMING_STEPS):
        test_pos, test_vel, test_acc, _ = velocity_verlet_step(
            test_pos, test_vel, test_acc, DT
        )
    t1 = time.perf_counter()
    
    ms_per_step = ((t1 - t0) / TIMING_STEPS) * 1000
    n_pairs = N * (N - 1) // 2
    timings.append(ms_per_step)
    print(f"{N:>6} | {ms_per_step:>16.4f} | {n_pairs:>10,}")

# --- PLOT ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: Raw timing
axes[0].plot(particle_counts, timings, 'o-', color='steelblue', lw=2, ms=8)
axes[0].set_xlabel("Number of Particles (N)", fontsize=12)
axes[0].set_ylabel("Wall Time per Step (ms)", fontsize=12)
axes[0].set_title("Raw Scaling: Time vs. N", fontsize=13)
axes[0].grid(True, alpha=0.4)

# Right: Log-log plot to reveal the power law
axes[1].loglog(particle_counts, timings, 'o-', color='crimson', lw=2, ms=8, label='Measured')

# Overlay a reference O(N^2) line
N_arr = np.array(particle_counts, dtype=float)
ref = timings[0] * (N_arr / particle_counts[0])**2
axes[1].loglog(particle_counts, ref, 'k--', lw=1.5, alpha=0.6, label=r'$O(N^2)$ reference')

axes[1].set_xlabel("Number of Particles (N)", fontsize=12)
axes[1].set_ylabel("Wall Time per Step (ms)", fontsize=12)
axes[1].set_title(r"Log-Log Scaling (slope ≈ 2 confirms $O(N^2)$)", fontsize=13)
axes[1].legend()
axes[1].grid(True, which='both', alpha=0.3)

plt.suptitle("compute_forces Scaling Analysis", fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.show()

print("\nConclusion: The log-log slope ≈ 2 confirms O(N²) scaling.")
print("Real codes use neighbor lists + cutoffs to reduce this to ~O(N).")
```

### What Does This Mean in Practice?

The table below uses your measured `time_per_step` to project realistic runtimes. Real MD codes (LAMMPS, GROMACS) use **neighbor lists** and parallel GPUs to achieve close to $O(N)$ scaling — but the same physics applies.

```python
# --- PROJECTION TABLE ---
# Uses the per-step time measured in Part A

print(f"\nProjection based on N={N_PARTICLES} atoms, {time_per_step_ms:.3f} ms/step\n")
print(f"{'Target simulation':25} | {'Steps needed':>14} | {'Estimated time':>16}")
print("-" * 60)

scenarios = [
    ("100 steps (this tutorial)",   100),
    ("1,000 steps",                 1_000),
    ("100,000 steps",               100_000),
    ("1 µs (1,000,000 steps)",      1_000_000),
    ("1 ms (1,000,000,000 steps)",  1_000_000_000),
]

for label, n_steps in scenarios:
    total_s = (time_per_step_ms / 1000) * n_steps
    if total_s < 60:
        time_str = f"{total_s:.2f} sec"
    elif total_s < 3600:
        time_str = f"{total_s/60:.1f} min"
    elif total_s < 86400:
        time_str = f"{total_s/3600:.1f} hours"
    else:
        time_str = f"{total_s/86400:.0f} days"
    print(f"{label:25} | {n_steps:>14,} | {time_str:>16}")

print("\nThis is why HPC clusters running thousands of CPU cores in")
print("parallel — like Anvil — are essential for real research!")
```
