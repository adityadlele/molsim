
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
