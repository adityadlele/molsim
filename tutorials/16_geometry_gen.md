# Tutorial 16: Geometry File Generation for LAMMPS

**Objective:** Build LAMMPS-ready geometry files for four classes of systems: metals, alloys, molecules, and polymers. Understand the LAMMPS data file format well enough to construct and troubleshoot input structures.

:::{note} Prerequisites
- Tutorial 4 (Python/Conda on Anvil)
- Basic ASE familiarity (Tutorial 13)
- The `molsimclass` conda environment
:::

---

## Setup

```bash
cd $SCRATCH
mkdir -p tutorial16
cd tutorial16
```

```python
import numpy as np
from ase import Atoms
from ase.build import bulk, molecule, make_supercell
from ase.io import write, read
import matplotlib.pyplot as plt

print("Imports successful.")
```

---

## Part 1: The LAMMPS Data File Format

Before generating files, it helps to understand what LAMMPS expects. A LAMMPS data file has a fixed structure:

```
# Comment line (optional)

N_atoms atoms
N_bonds bonds         (omit if no bonds)
N_angles angles       (omit if no angles)

N_atom_types atom types
N_bond_types bond types    (omit if no bonds)

xlo xhi xbounds
ylo yhi ybounds
zlo zhi zbounds

Masses

1  mass_type1
2  mass_type2

Atoms  # style_keyword

...one line per atom...

Bonds  (omit if no bonds)

...one line per bond...
```

### Atom Styles

The **atom style** determines what information each atom line contains. The most common styles are:

| Style | Fields per atom line | Use case |
|:------|:--------------------|:---------|
| `atomic` | id type x y z | Metals, noble gases |
| `charge` | id type charge x y z | Ionic materials, polarizable systems |
| `full` | id mol_id type charge x y z | Molecules with charges |
| `bond` | id mol_id type x y z | Polymers without charges |
| `molecular` | id mol_id type x y z | Alias for bond style |

The atom style must match the `atom_style` keyword in your LAMMPS input script. Mixing styles causes errors.

### Viewing a Data File

```python
# Read any LAMMPS data file back into ASE to verify it
atoms = read("output.lammps", format="lammps-data")
print(f"Read {len(atoms)} atoms")
print(f"Cell:\n{atoms.get_cell()}")
print(f"Atom types: {set(atoms.get_chemical_symbols())}")
```

---

## Part 2: Metals and Crystal Structures

### 2.1 The Three Common Crystal Structures

Most engineering metals adopt one of three crystal structures. The structure determines packing density, number of nearest neighbors, and slip systems — which directly controls ductility and strength.

**FCC (Face-Centered Cubic)**

Atoms sit at the corners and face centers of a cube. Each atom has **12 nearest neighbors** — the most densely packed arrangement possible for equal spheres.

- Lattice parameter: $a$ (one parameter fully defines the structure)
- Packing fraction: 74%
- Examples: Cu, Al, Au, Ag, Ni, Pb
- Slip system: {111}⟨110⟩ — 12 independent slip systems → **ductile**

**BCC (Body-Centered Cubic)**

Atoms at cube corners plus one in the center. Each atom has **8 nearest neighbors**.

- Packing fraction: 68%
- Examples: Fe (room T), W, Cr, Mo, V, Nb
- Slip system: {110}⟨111⟩ — fewer than FCC → harder, less ductile
- BCC metals often show a ductile-to-brittle transition at low temperature

**HCP (Hexagonal Close-Packed)**

Two alternating layers of hexagonally arranged atoms (ABABAB... stacking), compared to FCC's ABCABC stacking. Each atom has **12 nearest neighbors** — same as FCC.

- Two parameters: $a$ (in-plane spacing) and $c$ (layer spacing); ideal $c/a = \sqrt{8/3} \approx 1.633$
- Examples: Ti, Mg, Zn, Co, Zr
- Fewer active slip systems than FCC → generally less ductile than FCC metals

```python
# Visualize the three structures
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
structure_names = ['FCC (Cu)', 'BCC (Fe)', 'HCP (Ti)']

# Build each structure
fcc_cu = bulk('Cu', crystalstructure='fcc', a=3.615, cubic=True)
bcc_fe = bulk('Fe', crystalstructure='bcc', a=2.870, cubic=True)
hcp_ti = bulk('Ti', crystalstructure='hcp', a=2.951, c=4.684)

structures = [fcc_cu, bcc_fe, hcp_ti]

for ax, atoms, name in zip(axes, structures, structure_names):
    pos = atoms.get_positions()
    ax.scatter(pos[:, 0], pos[:, 1], s=300, c='steelblue',
               edgecolors='black', linewidths=1.5, zorder=5)
    cell = atoms.get_cell()
    # Draw cell outline (projected onto xy plane)
    corners = [(0,0), (cell[0,0], cell[0,1]),
               (cell[0,0]+cell[1,0], cell[0,1]+cell[1,1]),
               (cell[1,0], cell[1,1]), (0,0)]
    xs, ys = zip(*corners)
    ax.plot(xs, ys, 'k-', lw=1.5, alpha=0.5)
    ax.set_title(f"{name}\n{len(atoms)} atoms/cell", fontsize=11)
    ax.set_aspect('equal')
    ax.axis('off')

plt.suptitle("Unit Cells: FCC, BCC, HCP (xy projection)", fontsize=13)
plt.tight_layout()
plt.savefig("crystal_structures.png", dpi=150)
plt.show()
```

### 2.2 Generating Metal Data Files

```python
# ── FCC Copper ────────────────────────────────────────────────────
cu = bulk('Cu', crystalstructure='fcc', a=3.615, cubic=True)
cu_super = make_supercell(cu, [[5, 0, 0], [0, 5, 0], [0, 0, 5]])
write("cu_fcc.lammps", cu_super, format="lammps-data",
      atom_style="atomic")

print(f"Cu FCC supercell: {len(cu_super)} atoms")
print(f"Cell dimensions: {cu_super.get_cell().diagonal().round(3)} Å")

# ── BCC Iron ──────────────────────────────────────────────────────
fe = bulk('Fe', crystalstructure='bcc', a=2.870, cubic=True)
fe_super = make_supercell(fe, [[5, 0, 0], [0, 5, 0], [0, 0, 5]])
write("fe_bcc.lammps", fe_super, format="lammps-data",
      atom_style="atomic")

print(f"Fe BCC supercell: {len(fe_super)} atoms")

# ── HCP Titanium ──────────────────────────────────────────────────
ti = bulk('Ti', crystalstructure='hcp', a=2.951, c=4.684)
# HCP supercell: repeat in all three directions
ti_super = ti.repeat([6, 6, 4])
write("ti_hcp.lammps", ti_super, format="lammps-data",
      atom_style="atomic")

print(f"Ti HCP supercell: {len(ti_super)} atoms")
```

### 2.3 Verifying the Output

```python
# Read back and check
for fname, symbol in [("cu_fcc.lammps", "Cu"),
                       ("fe_bcc.lammps", "Fe"),
                       ("ti_hcp.lammps", "Ti")]:
    atoms = read(fname, format="lammps-data",
                 style="atomic", units="metal")
    symbols = atoms.get_chemical_symbols()
    cell    = atoms.get_cell().diagonal()
    print(f"{fname}: {len(atoms)} atoms, "
          f"cell = {cell.round(2)} Å, "
          f"elements = {set(symbols)}")
```

### 2.4 Corresponding LAMMPS Input

```lammps
# lammps_cu.in — use with cu_fcc.lammps
units       metal
atom_style  atomic
boundary    p p p

read_data   cu_fcc.lammps

pair_style  eam
pair_coeff  * * /anvil/projects/x-chm250117/potentials/Cu_u3.eam

thermo      100
run         0    # Just verify structure loads; add dynamics as needed
```

---

## Part 3: Alloys

### 3.1 Random Substitutional Alloy (Cu-Au)

A solid solution where Au atoms randomly replace some Cu atoms on the FCC lattice. The composition is set by the substitution fraction.

```python
import random

def random_alloy(base_symbol, sub_symbol, base_a, composition,
                 supercell_size=5, seed=42):
    """
    Create a random substitutional binary alloy.

    Parameters
    ----------
    base_symbol  : str   — majority element symbol
    sub_symbol   : str   — minority element symbol
    base_a       : float — lattice constant of the base metal (Å)
    composition  : float — fraction of sub_symbol atoms (0–1)
    supercell_size: int  — supercell repeat in each direction
    seed         : int   — random seed for reproducibility

    Returns
    -------
    atoms : ASE Atoms object
    """
    random.seed(seed)

    # Start with pure base metal supercell
    unit = bulk(base_symbol, crystalstructure='fcc', a=base_a, cubic=True)
    atoms = make_supercell(unit, [[supercell_size, 0, 0],
                                   [0, supercell_size, 0],
                                   [0, 0, supercell_size]])

    # Choose which atoms to substitute
    n_total = len(atoms)
    n_sub   = int(round(n_total * composition))
    sub_indices = random.sample(range(n_total), n_sub)

    # Replace atom symbols at those indices
    symbols = atoms.get_chemical_symbols()
    for i in sub_indices:
        symbols[i] = sub_symbol
    atoms.set_chemical_symbols(symbols)

    return atoms


# Cu-10at%Au random alloy
cu_au = random_alloy('Cu', 'Au', base_a=3.615, composition=0.10)

n_cu = cu_au.get_chemical_symbols().count('Cu')
n_au = cu_au.get_chemical_symbols().count('Au')
print(f"Cu-Au alloy: {len(cu_au)} atoms total")
print(f"  Cu: {n_cu}  ({100*n_cu/len(cu_au):.1f} at%)")
print(f"  Au: {n_au}  ({100*n_au/len(cu_au):.1f} at%)")

write("cu_au_alloy.lammps", cu_au, format="lammps-data",
      atom_style="atomic")
print("Written: cu_au_alloy.lammps")
```

### 3.2 Understanding Multi-Type Data Files

When the alloy has two elements, the data file has two atom types. The LAMMPS input must assign the correct potential to each type:

```lammps
# lammps_cu_au.in
units       metal
atom_style  atomic
boundary    p p p

read_data   cu_au_alloy.lammps

# EAM/alloy potential — must list ALL types in order
pair_style  eam/alloy
pair_coeff  * * /path/to/CuAu.eam.alloy Cu Au
#                                        type1 type2
```

The order of element names in `pair_coeff` must match the atom type numbering in the data file. Type 1 in the data file gets the first element name, type 2 gets the second.

### 3.3 Ordered Intermetallic: L1₂ Structure (Cu₃Au)

Below ~390°C, the Cu-25%Au alloy orders into the L1₂ structure: Au atoms occupy cube corners, Cu atoms occupy face centers. This is a classic ordered intermetallic.

```python
# L1₂ Cu₃Au: manually set atom positions
a = 3.74   # Å (experimental for ordered Cu₃Au)

# 4-atom unit cell: 1 Au at corner, 3 Cu at face centers
positions = [
    (0.0, 0.0, 0.0),                  # Au at corner
    (0.5, 0.5, 0.0),                  # Cu at face
    (0.5, 0.0, 0.5),                  # Cu at face
    (0.0, 0.5, 0.5),                  # Cu at face
]
symbols = ['Au', 'Cu', 'Cu', 'Cu']

cu3au_unit = Atoms(
    symbols   = symbols,
    positions = [(np.array(p) * a) for p in positions],
    cell      = [a, a, a],
    pbc       = True
)

cu3au_super = make_supercell(cu3au_unit,
                              [[4, 0, 0], [0, 4, 0], [0, 0, 4]])

n_au = cu3au_super.get_chemical_symbols().count('Au')
n_cu = cu3au_super.get_chemical_symbols().count('Cu')
print(f"L1₂ Cu₃Au: {len(cu3au_super)} atoms")
print(f"  Cu: {n_cu}  Au: {n_au}  ratio: {n_cu/n_au:.0f}:1")

write("cu3au_ordered.lammps", cu3au_super, format="lammps-data",
      atom_style="atomic")
print("Written: cu3au_ordered.lammps")
```

---

## Part 4: Molecules

### 4.1 Writing LAMMPS Data Files with Topology

Molecules require **bonds, angles, and dihedrals** in the data file — these define the molecular topology that LAMMPS uses to apply harmonic bond and angle potentials. ASE's `lammps-data` writer handles atomic systems but not full molecular topology. We write a Python function for this.

### 4.2 Water Molecule (SPC/E Model)

The SPC/E water model places a negative charge on oxygen and positive charges on hydrogen. A box of water molecules starts with a single optimized geometry that is then replicated.

```python
def write_water_lammps(n_molecules, box_length, filename,
                        seed=42):
    """
    Write a box of SPC/E water molecules as a LAMMPS data file.

    SPC/E charges: O = -0.8476 e,  H = +0.4238 e
    Geometry:      O-H bond = 1.0 Å,  H-O-H angle = 109.47°

    Parameters
    ----------
    n_molecules : int   — number of water molecules
    box_length  : float — cubic box side length (Å)
    filename    : str   — output file path
    seed        : int   — random seed for placement
    """
    rng = np.random.default_rng(seed)

    # SPC/E geometry
    r_OH    = 1.0      # Å
    theta   = np.radians(109.47)

    # Generate positions: O at random locations, H placed relative to O
    atoms_out  = []   # (type, charge, x, y, z)
    bonds_out  = []   # (bond_type, i, j)  1-indexed
    angles_out = []   # (angle_type, i, j, k)  j is vertex

    for m in range(n_molecules):
        O_pos = rng.uniform(1.5, box_length - 1.5, size=3)

        # Random orientation for H atoms
        phi   = rng.uniform(0, 2 * np.pi)
        H1_pos = O_pos + r_OH * np.array([np.sin(theta/2) * np.cos(phi),
                                           np.sin(theta/2) * np.sin(phi),
                                           np.cos(theta/2)])
        H2_pos = O_pos + r_OH * np.array([-np.sin(theta/2) * np.cos(phi),
                                           -np.sin(theta/2) * np.sin(phi),
                                            np.cos(theta/2)])

        base = 3 * m + 1    # 1-indexed atom counter
        atoms_out.append((1, -0.8476, *O_pos))   # O
        atoms_out.append((2,  0.4238, *H1_pos))  # H
        atoms_out.append((2,  0.4238, *H2_pos))  # H

        bonds_out.append((1, base,   base+1))     # O-H1
        bonds_out.append((1, base,   base+2))     # O-H2
        angles_out.append((1, base+1, base, base+2))  # H-O-H

    n_atoms  = len(atoms_out)
    n_bonds  = len(bonds_out)
    n_angles = len(angles_out)

    with open(filename, 'w') as f:
        f.write(f"# SPC/E water  {n_molecules} molecules\n\n")
        f.write(f"{n_atoms} atoms\n")
        f.write(f"{n_bonds} bonds\n")
        f.write(f"{n_angles} angles\n\n")
        f.write("2 atom types\n")
        f.write("1 bond types\n")
        f.write("1 angle types\n\n")
        f.write(f"0.0 {box_length:.4f} xlo xhi\n")
        f.write(f"0.0 {box_length:.4f} ylo yhi\n")
        f.write(f"0.0 {box_length:.4f} zlo zhi\n\n")
        f.write("Masses\n\n")
        f.write("1 15.9994  # O\n")
        f.write("2  1.0080  # H\n\n")
        f.write("Atoms  # full: id mol_id type charge x y z\n\n")
        for i, (atype, charge, x, y, z) in enumerate(atoms_out, start=1):
            mol_id = (i - 1) // 3 + 1
            f.write(f"{i} {mol_id} {atype} {charge:.4f} "
                    f"{x:.6f} {y:.6f} {z:.6f}\n")
        f.write("\nBonds\n\n")
        for i, (btype, a1, a2) in enumerate(bonds_out, start=1):
            f.write(f"{i} {btype} {a1} {a2}\n")
        f.write("\nAngles\n\n")
        for i, (atype, a1, a2, a3) in enumerate(angles_out, start=1):
            f.write(f"{i} {atype} {a1} {a2} {a3}\n")

    print(f"Written: {filename}")
    print(f"  {n_molecules} molecules, {n_atoms} atoms, "
          f"{n_bonds} bonds, {n_angles} angles")


# Generate a small water box
write_water_lammps(n_molecules=100, box_length=14.7, filename="water_100.lammps")
```

The corresponding LAMMPS input:

```lammps
# water_spc_e.in
units       real
atom_style  full
boundary    p p p

read_data   water_100.lammps

pair_style  lj/cut/coul/long 10.0
pair_coeff  1 1 0.1553 3.1660   # O-O LJ
pair_coeff  2 2 0.0    0.0      # H-H (no LJ)

bond_style  harmonic
bond_coeff  1 1000.0 1.0        # stiff: O-H rigid approximation

angle_style harmonic
angle_coeff 1 100.0 109.47      # H-O-H angle

kspace_style pppm 1.0e-4        # long-range electrostatics
```

---

## Part 5: Polymers

### 5.1 Polyethylene Chain Builder

A polymer chain is a sequence of monomer units connected by covalent bonds. For polyethylene (–CH₂–)ₙ, we build a chain in the united-atom model: each CH₂ group is a single "bead" with the combined mass of C + 2H.

```python
def build_polyethylene_lammps(n_chains, n_monomers_per_chain,
                               box_length, filename,
                               bond_length=1.54, seed=42):
    """
    Build a polyethylene system in the united-atom (UA) model.

    Each CH₂ or CH₃ group is one bead (atom type).
    Atom types:
      1 = CH₃ (chain end, mass 15.035 amu)
      2 = CH₂ (chain middle, mass 14.027 amu)

    Bond: CH₂–CH₂ at 1.54 Å
    Angle: CH₂–CH₂–CH₂ at 114°

    Parameters
    ----------
    n_chains           : int   — number of polymer chains
    n_monomers_per_chain: int  — degree of polymerization (chain length)
    box_length         : float — cubic box side length (Å)
    filename           : str   — output file path
    bond_length        : float — C-C bond length (Å)
    seed               : int   — random seed
    """
    rng = np.random.default_rng(seed)
    theta_eq = np.radians(114.0)  # equilibrium C-C-C angle

    atoms_list  = []   # (mol_id, type, x, y, z)
    bonds_list  = []   # (bond_type, i, j)
    angles_list = []   # (angle_type, i, j, k)

    atom_counter = 0

    for chain in range(n_chains):
        mol_id = chain + 1

        # Place first bead randomly
        pos = [rng.uniform(1.0, box_length - 1.0, size=3)]

        # Grow chain bead by bead using random walk with fixed bond angle
        direction = rng.uniform(-1, 1, size=3)
        direction /= np.linalg.norm(direction)

        for i in range(1, n_monomers_per_chain):
            # Perturb direction slightly (restricted random walk)
            perp = rng.uniform(-1, 1, size=3)
            perp -= np.dot(perp, direction) * direction
            if np.linalg.norm(perp) > 1e-10:
                perp /= np.linalg.norm(perp)

            # New direction at the prescribed bond angle
            new_dir = (np.cos(np.pi - theta_eq) * direction +
                       np.sin(np.pi - theta_eq) * perp)
            new_dir /= np.linalg.norm(new_dir)
            direction = new_dir

            new_pos = pos[-1] + bond_length * direction
            # Wrap into box (simple approach)
            new_pos = new_pos % box_length
            pos.append(new_pos)

        # Add atoms: type 1 for end beads, type 2 for middle beads
        chain_start = atom_counter + 1
        for i, p in enumerate(pos):
            atom_counter += 1
            atype = 1 if (i == 0 or i == len(pos) - 1) else 2
            atoms_list.append((mol_id, atype, p[0], p[1], p[2]))

        # Add bonds within this chain
        for i in range(n_monomers_per_chain - 1):
            bonds_list.append((1, chain_start + i, chain_start + i + 1))

        # Add angles within this chain
        for i in range(n_monomers_per_chain - 2):
            angles_list.append((1, chain_start + i,
                                    chain_start + i + 1,
                                    chain_start + i + 2))

    n_atoms  = len(atoms_list)
    n_bonds  = len(bonds_list)
    n_angles = len(angles_list)

    with open(filename, 'w') as f:
        f.write(f"# Polyethylene UA: {n_chains} chains × "
                f"{n_monomers_per_chain} monomers\n\n")
        f.write(f"{n_atoms}  atoms\n")
        f.write(f"{n_bonds}  bonds\n")
        f.write(f"{n_angles}  angles\n\n")
        f.write("2 atom types\n")
        f.write("1 bond types\n")
        f.write("1 angle types\n\n")
        f.write(f"0.0 {box_length:.4f} xlo xhi\n")
        f.write(f"0.0 {box_length:.4f} ylo yhi\n")
        f.write(f"0.0 {box_length:.4f} zlo zhi\n\n")
        f.write("Masses\n\n")
        f.write("1 15.035  # CH3 end bead\n")
        f.write("2 14.027  # CH2 middle bead\n\n")
        f.write("Atoms  # bond: id mol_id type x y z\n\n")
        for i, (mol_id, atype, x, y, z) in enumerate(atoms_list, start=1):
            f.write(f"{i} {mol_id} {atype} {x:.6f} {y:.6f} {z:.6f}\n")
        f.write("\nBonds\n\n")
        for i, (btype, a1, a2) in enumerate(bonds_list, start=1):
            f.write(f"{i} {btype} {a1} {a2}\n")
        f.write("\nAngles\n\n")
        for i, (atype, a1, a2, a3) in enumerate(angles_list, start=1):
            f.write(f"{i} {atype} {a1} {a2} {a3}\n")

    print(f"Written: {filename}")
    print(f"  {n_chains} chains × {n_monomers_per_chain} monomers = "
          f"{n_atoms} atoms total")
    print(f"  {n_bonds} bonds,  {n_angles} angles")


# Generate polyethylene: 10 chains of 50 monomers each
build_polyethylene_lammps(
    n_chains=10,
    n_monomers_per_chain=50,
    box_length=40.0,
    filename="polyethylene_10x50.lammps"
)
```

The corresponding LAMMPS input (TraPPE united-atom force field):

```lammps
# polyethylene_ua.in
units       real
atom_style  bond
boundary    p p p

read_data   polyethylene_10x50.lammps

pair_style  lj/cut 14.0
pair_coeff  1 1 0.1947 3.75   # CH3-CH3 TraPPE
pair_coeff  1 2 0.1730 3.83   # CH3-CH2
pair_coeff  2 2 0.0914 3.95   # CH2-CH2

bond_style  harmonic
bond_coeff  1 350.0 1.54      # C-C stretch (stiff)

angle_style harmonic
angle_coeff 1  62.5 114.0     # C-C-C TraPPE angle

special_bonds lj 0.0 0.0 1.0  # exclude 1-2 and 1-3 non-bonded pairs

thermo      1000
run         0
```

---

## Part 6: Common Pitfalls and Debugging

**Atom type mismatch:** The number of atom types declared in the header must match the number of `pair_coeff` lines in the LAMMPS input. A 2-type data file needs coefficients for types 1-1, 1-2, and 2-2 (if using a pair style).

**Unit consistency:** ASE writes `lammps-data` files in whatever units the atoms object uses internally (typically Å). Your LAMMPS input must use `units metal` (for metals with Å/eV) or `units real` (for molecules with Å/kcal/mol). Mixing unit systems is a common and silent source of wrong results.

**Box too small:** For a valid periodic simulation, each box dimension must be at least $2 \times r_{\text{cutoff}}$. For a 10 Å cutoff, each dimension must be at least 20 Å. Smaller boxes cause atoms to interact with their own periodic images.

**Overlapping atoms:** Random placement of molecules can generate atom-atom distances smaller than the repulsive core of the potential, causing forces to become enormous and the simulation to explode in the first few steps. Always verify minimum distances after building a structure:

```python
from ase.geometry import get_distances

atoms = read("water_100.lammps", format="lammps-data",
             style="full", units="real")
pos   = atoms.get_positions()
cell  = atoms.get_cell()

# Check all pairwise distances (small system only)
dists, _ = get_distances(pos, cell=cell, pbc=True)
np.fill_diagonal(dists, np.inf)
print(f"Minimum interatomic distance: {dists.min():.3f} Å")
# Should be > ~1.5 Å for any reasonable initial configuration
```

**Polymer chain crossing box boundaries:** The simple modulo wrapping in the polymer builder above can create bonds that span most of the box. For production use, tools like Packmol or Moltemplate handle chain placement much more carefully.

---

## Summary

| System | ASE builder | Atom style | Key file sections |
|:-------|:-----------|:-----------|:-----------------|
| Metal (FCC/BCC) | `bulk()` + `make_supercell()` | `atomic` | Masses, Atoms |
| Metal (HCP) | `bulk(..., hcp)` + `repeat()` | `atomic` | Masses, Atoms |
| Random alloy | `bulk()` + random symbol swap | `atomic` | Masses, Atoms (2 types) |
| Ordered intermetallic | Manual `Atoms()` | `atomic` | Masses, Atoms (2 types) |
| Molecule | Custom writer | `full` | Masses, Atoms, Bonds, Angles |
| Polymer | Custom chain builder | `bond` | Masses, Atoms, Bonds, Angles |

The two custom writers in Parts 4 and 5 are starting points. For production molecular and polymer simulations, consider **Moltemplate** (for complex force fields) or **Packmol** (for dense liquid/amorphous packing) — both generate LAMMPS-ready data files and handle edge cases that simple Python builders miss.

---

## Further Reading

- LAMMPS data file format: [https://docs.lammps.org/read_data.html](https://docs.lammps.org/read_data.html)
- ASE I/O formats: [https://wiki.fysik.dtu.dk/ase/ase/io/io.html](https://wiki.fysik.dtu.dk/ase/ase/io/io.html)
- Moltemplate: [https://www.moltemplate.org](https://www.moltemplate.org)
- Packmol: [http://m3g.iqm.unicamp.br/packmol](http://m3g.iqm.unicamp.br/packmol)
- TraPPE force field for alkanes: [http://trappe.oit.umn.edu](http://trappe.oit.umn.edu)
