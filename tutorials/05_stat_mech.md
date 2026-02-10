
# Tutorial 5: Analyzing Argon MD Simulations

## Introduction
In this tutorial, we will process the raw output of a Molecular Dynamics (MD) simulation of Argon gas.

**The Goal:**
We want to verify if our microscopic simulation correctly reproduces the macroscopic properties we intended. specifically:
1.  Does the system maintain the target temperature ($300 \text{ K}$)?
2.  Do the atoms follow the correct velocity distribution (Maxwell-Boltzmann)?

### Prerequisites
* Ensure you have the file `velocities.txt` in your working directory.
* This file is a LAMMPS dump containing the velocities ($v_x, v_y, v_z$) of Argon atoms over several timesteps.

---

## 1. Setup and Unit Conversions

MD simulations often use "real" units (Angstroms, femtoseconds). However, to use the Boltzmann constant in standard SI units ($J/K$), we must convert our velocities to meters per second ($m/s$).

**Copy and run the following code to set up our constants:**

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import maxwell

# --- 1. Constants & Unit Conversions ---
KB = 1.380649e-23        # Boltzmann constant (J/K)
NA = 6.02214076e23       # Avogadro's number (mol^-1)
MASS_AR_G_MOL = 39.95    # Mass of Argon (g/mol)

# Convert Mass to kg per atom
# (g/mol -> kg/mol -> kg/atom)
mass_kg = (MASS_AR_G_MOL * 1e-3) / NA

# Conversion factor for Velocity
# LAMMPS 'real' units: Angstroms/femtosecond
# 1 Angstrom = 1e-10 m
# 1 fs = 1e-15 s
# Factor = 1e-10 / 1e-15 = 1e5
VELOCITY_CONV_FACTOR = 1e5 

```

---

## 2. Parsing the LAMMPS Dump File

The simulation data is stored in a text file where each "snapshot" of time is listed sequentially. We need a function to read this file and extract the velocity vectors () for every atom at every timestep.

**Copy this function to parse the data:**

```python
def parse_lammps_dump(filename):
    """
    Parses a LAMMPS dump file and returns a list of time steps 
    and a list of velocity arrays.
    """
    timesteps = []
    velocities = []
    
    current_velocities = []
    reading_atoms = False
    
    with open(filename, 'r') as f:
        for line in f:
            if "ITEM: TIMESTEP" in line:
                if current_velocities:
                    velocities.append(np.array(current_velocities))
                    current_velocities = []
                # Read the next line which is the timestep number
                timesteps.append(int(next(f)))
                reading_atoms = False
            
            elif "ITEM: ATOMS" in line:
                reading_atoms = True
                continue
                
            elif reading_atoms:
                # Line format: id type vx vy vz
                parts = line.split()
                # We only need vx, vy, vz (indices 2, 3, 4)
                vx = float(parts[2])
                vy = float(parts[3])
                vz = float(parts[4])
                current_velocities.append([vx, vy, vz])
                
        # Append the last frame
        if current_velocities:
            velocities.append(np.array(current_velocities))
            
    return np.array(timesteps), velocities

# --- Load the Data ---
filename = 'velocities.txt'
print(f"Reading {filename}...")
steps, velocity_data = parse_lammps_dump(filename)
print(f"Successfully processed {len(steps)} frames.")

```

---

## 3. Calculating Temperature (The Physics)

Now we apply the **Equipartition Theorem**. We learned that temperature is directly related to the average Kinetic Energy () of the system:

$$ \langle KE \rangle = \frac{3}{2} k_B T \quad \Rightarrow \quad T = \frac{2 \langle KE \rangle}{3 k_B} $$

We will calculate the instantaneous temperature for every frame in our simulation.

```python
temperatures = []

for frame_vel in velocity_data:
    # frame_vel is an array of shape (N_atoms, 3)
    
    # 1. Convert velocity to m/s
    v_ms = frame_vel * VELOCITY_CONV_FACTOR
    
    # 2. Calculate squared speed (v^2 = vx^2 + vy^2 + vz^2)
    v_sq = np.sum(v_ms**2, axis=1)
    
    # 3. Calculate Kinetic Energy per atom (Joules)
    # KE = 1/2 * m * v^2
    ke_per_atom = 0.5 * mass_kg * v_sq
    
    # 4. Calculate Average Kinetic Energy for the system
    avg_ke = np.mean(ke_per_atom)
    
    # 5. Calculate Temperature from Equipartition Theorem
    temp = (2 * avg_ke) / (3 * KB)
    
    temperatures.append(temp)

# Calculate the Ensemble Average (mean over all time)
ensemble_avg_temp = np.mean(temperatures)

print("-" * 30)
print(f"Target Temperature: 300.00 K")
print(f"Calculated Temperature: {ensemble_avg_temp:.2f} K")
print("-" * 30)

```

---

## 4. Visualizing Temperature Fluctuations

In an **NVT Ensemble** (Constant Number, Volume, Temperature), the temperature is not perfectly fixed. It fluctuates around a mean value because the thermostat adds and removes energy to keep it stable.

Let's plot these fluctuations to see if our simulation has "equilibrated" (settled down).

```python
plt.figure(figsize=(10, 6))
plt.plot(steps, temperatures, marker='o', linestyle='-', color='b', label='Instantaneous T')
plt.axhline(y=300, color='r', linestyle='--', label='Thermostat Target (300K)')
plt.axhline(y=ensemble_avg_temp, color='g', linestyle='-', label=f'Calculated Mean ({ensemble_avg_temp:.1f} K)')

plt.xlabel('Time Step')
plt.ylabel('Temperature (K)')
plt.title('Temperature Fluctuation in NVT Ensemble (Argon)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()

```

---

## 5. Validating the Velocity Distribution

Finally, we perform the ultimate check. Do our atoms follow the laws of statistical mechanics?
We will plot a histogram of all atom speeds and compare it to the theoretical **Maxwell-Boltzmann** curve for Argon at 300K.

```python
# --- 1. Collect all speeds from all frames ---
all_speeds = []

for frame_vel in velocity_data:
    # Convert vectors to speeds (magnitude)
    # v = sqrt(vx^2 + vy^2 + vz^2)
    speeds_si = np.linalg.norm(frame_vel, axis=1) * VELOCITY_CONV_FACTOR
    all_speeds.extend(speeds_si)

all_speeds = np.array(all_speeds)

# --- 2. Define Theoretical Maxwell-Boltzmann Curve ---
T_target = 300.0
v_theory = np.linspace(0, 1000, 100)

# The Maxwell-Boltzmann probability density function formula
# f(v) = 4*pi * (m / (2*pi*k*T))^(3/2) * v^2 * exp(-m*v^2 / (2*k*T))
prefactor = 4 * np.pi * (mass_kg / (2 * np.pi * KB * T_target))**(1.5)
mb_curve = prefactor * (v_theory**2) * np.exp(-mass_kg * (v_theory**2) / (2 * KB * T_target))

# --- 3. Plotting ---
plt.figure(figsize=(10, 6))

# Histogram of simulation data
# 'density=True' normalizes the histogram so area under curve = 1
plt.hist(all_speeds, bins=50, density=True, alpha=0.6, color='skyblue', edgecolor='black', label='Simulation Data')

# Theoretical Curve
plt.plot(v_theory, mb_curve, 'r-', linewidth=2, label=f'Maxwell-Boltzmann ({T_target} K)')

plt.xlabel('Speed (m/s)')
plt.ylabel('Probability Density')
plt.title(f'Velocity Distribution of Argon at {T_target} K')
plt.legend()
plt.grid(True, alpha=0.3)
plt.xlim(0, 1000)

plt.show()

```

```

```