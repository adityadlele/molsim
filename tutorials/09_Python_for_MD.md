
# Tutorial 9: Analyzing LAMMPS Output with Python and Jupyter

Now that you have successfully run your first Molecular Dynamics simulation on Anvil, it is time to make sense of the results. LAMMPS produces text files containing massive amounts of raw data. To extract meaningful physics, we will use Python within a Jupyter Notebook to parse, visualize, and analyze this data.

---

## 1. Understanding the LAMMPS Log File

When LAMMPS runs, it generates a file typically named `log.lammps`. This file contains a record of everything that happened during your job, including:
1. The script commands it read.
2. System initialization details (number of atoms, box size).
3. The **thermodynamic output** (the data we want!).
4. Performance and timing statistics at the end.



The thermodynamic output block looks something like a spreadsheet printed in plain text. It starts with a header row matching the `thermo_style` you defined in your script:

```text
Step Temp E_pair E_mol TotEng Press 
       0          300   -1.5032014            0    11.144988    63.774488 
     100    274.69837   -1.2268808            0    10.368146   -60.627771 
     200    222.13327   -1.2384931            0     8.139413   -36.938362 
...
Loop time of 0.1234 on 4 procs for 10000 steps

```

Our goal is to extract this block of numbers into a Python Pandas DataFrame.

---

## 2. Parsing the Data with Python

Instead of copying and pasting the numbers manually (which is impossible for large simulations!), we can write a small Python function to automatically find the data block and load it into a DataFrame.

Open a new Jupyter Notebook and run the following code cell to import the necessary libraries and define our parser function:

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def read_lammps_log(filepath):
    """
    Reads a LAMMPS log file and extracts the thermodynamic data into a Pandas DataFrame.
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()
        
    start_idx = -1
    end_idx = -1
    
    # Scan through the file to find where the data starts and ends
    for i, line in enumerate(lines):
        # Look for the header row (adjust this if your thermo_style changes!)
        if line.startswith("Step Temp"): 
            start_idx = i
        # Look for the end of the run
        elif start_idx != -1 and line.startswith("Loop time"):
            end_idx = i
            break
            
    if start_idx == -1 or end_idx == -1:
        print("Error: Could not find the thermodynamic data block.")
        return None
        
    # Extract just the rows containing the data
    data_lines = lines[start_idx:end_idx]
    
    # Split each row into columns
    parsed_data = [row.split() for row in data_lines]
    
    # Convert to a Pandas DataFrame (first row is headers, the rest is data)
    df = pd.DataFrame(parsed_data[1:], columns=parsed_data[0], dtype=float)
    
    return df

# Load the data! (Make sure log.lammps is in the same folder as your notebook)
df = read_lammps_log("log.lammps")

# Print the first 5 rows to verify it worked
print(df.head())

```

---

## 3. Visualizing Equilibration vs. Production

In Molecular Dynamics, the initial state of the system is often artificial (e.g., atoms placed on a perfect grid, or given exact velocities). The system must run for a period of time to "relax" into a natural state. This is called the **equilibration phase**. We only calculate macroscopic properties from the **production phase** that follows.

Let's plot the Temperature and Total Energy to visually identify these phases.

### Plotting Temperature

```python
plt.figure(figsize=(8, 5))
plt.plot(df['Step'], df['Temp'], color='red', label='Instantaneous Temp')

plt.axvline(x=2000, color='black', linestyle='--', label='End of Equilibration')

plt.xlabel('Time Step')
plt.ylabel('Temperature (K)')
plt.title('System Temperature over Time')
plt.legend()
plt.grid(True)
plt.show()

```

Notice how the temperature might spike or drop dramatically at the beginning before settling into stable fluctuations around your target value.

### Plotting Energy (Checking Physics)

If you ran an **NVE** simulation, the Total Energy should be a flat horizontal line, even though the Kinetic and Potential energies are fluctuating.

```python
plt.figure(figsize=(8, 5))
plt.plot(df['Step'], df['pe'], label='Potential Energy (PE)')
plt.plot(df['Step'], df['ke'], label='Kinetic Energy (KE)')
plt.plot(df['Step'], df['etotal'], label='Total Energy (PE + KE)', color='black', linewidth=2)

plt.xlabel('Time Step')
plt.ylabel('Energy (kcal/mol)')
plt.title('Energy Conservation in NVE Ensemble')
plt.legend()
plt.grid(True)
plt.show()

```

---

## 4. Calculating Macroscopic Properties

Once you have identified where your production phase begins, you can slice your DataFrame to calculate the average equilibrium properties. Let's say we look at our plots and decide the system is fully equilibrated after step 2000.

```python
# Create a new DataFrame containing only the production data
# (Keep only rows where the Step is greater than 2000)
production_df = df[df['Step'] > 2000]

# Calculate Mean and Standard Deviation for Pressure
mean_pressure = production_df['Press'].mean()
std_pressure = production_df['Press'].std()

# Calculate Mean for Temperature
mean_temp = production_df['Temp'].mean()

print(f"--- Equilibrium Properties ---")
print(f"Average Temperature: {mean_temp:.2f} K")
print(f"Average Pressure:    {mean_pressure:.2f} +/- {std_pressure:.2f} atm")

```

