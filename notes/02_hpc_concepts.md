# High Performance Computing (HPC) Architecture

Before logging in, it is critical to understand that a Supercomputer (HPC) is not just a "fast laptop." It is a fundamentally different architecture designed for parallel tasks.

In this course, we will use **[Purdue Anvil](https://www.rcac.purdue.edu/compute/anvil)** as our primary resource.

## What is an HPC Cluster?
Imagine you need to move a house. You have two options:
1.  **Vertical Scaling:** Buy a giant, custom-built truck the size of a house. (Expensive, rare).
2.  **Horizontal Scaling:** Rent 50 standard pickup trucks and coordinate them. (Cost-effective, scalable).

An HPC Cluster is option #2. It is a collection of thousands of individual computers (called **Nodes**) connected by a high-speed network.

### PC vs. HPC Comparison

| Feature | Your Laptop | HPC Node (Anvil) |
| :--- | :--- | :--- |
| **CPU** | Intel i7 / M1 (4-12 Cores) | [AMD EPYC "Milan"](https://www.rcac.purdue.edu/knowledge/anvil/specs) (128 Cores) |
| **RAM** | 16 GB - 32 GB | 256 GB - 1000 GB |
| **OS** | Windows / macOS | Linux (RHEL/CentOS) |
| **Usage** | Interactive (Mouse/GUI) | Batch (Scripts/Command Line) |

:::{seealso}
For a deep dive into the hardware, read the [Official Anvil System Specifications](https://www.rcac.purdue.edu/knowledge/anvil/specs).
:::

## The "Secret Sauce": The Interconnect
Why can't we just connect 100 laptops with Wi-Fi?

**Latency.**
In a molecular dynamics simulation, atoms on Node 1 interact with atoms on Node 2. They must exchange data every single time step ($10^{-15}$ seconds). If the network is slow, the CPUs sit idle waiting for data.

Anvil uses **NVIDIA HDR InfiniBand**.
* **Latency:** < 1 microsecond (vs 1000+ $\mu$s for Ethernet).
* **Bandwidth:** 100 Gigabits/sec.

## Storage Architecture: Where do I put my files?
This is the #1 reason new users lose data or crash simulations. HPC storage is divided into tiers with different rules.

Read the official [Anvil Storage Guide](https://www.rcac.purdue.edu/knowledge/anvil/storage) for exact quotas.

### 1. Home Directory (`$HOME`)
* **Location:** `/home/x-alele`
* **Purpose:** Configuration files (`.bashrc`), source code, small scripts.
* **Pros:** Backed up, safe.
* **Cons:** **Tiny** (usually 10-25 GB quota). **Slow** for I/O.
* **Warning:** *Never run a simulation here.* You will hit the quota instantly.

### 2. Scratch Directory (`$SCRATCH`)
* **Location:** `/anvil/scratch/x-alele`
* **Purpose:** Running simulations, storing trajectory files (`.dcd`, `.xtc`), large datasets.
* **Pros:** **Huge** (100 TB+), extremely fast parallel file system (Lustre).
* **Cons:** **Not Backed Up.** Files are **purged** (deleted) automatically if not touched for 90 days.

### 3. Project Directory (`$PROJECT`)
* **Location:** `/anvil/projects/class_allocation`
* **Purpose:** Shared data for the class, common datasets, or group collaboration.
* **Pros:** Persistent (no purge), shared access.
* **Cons:** Limited space compared to Scratch.

:::{danger} The Golden Rule of Storage
**Compute on Scratch, Backup to Home/Cloud.**
1.  Run your simulation in `$SCRATCH`.
2.  Analyze the data.
3.  Move the final results/graphs to Projects or `$HOME`.
4.  The system administrator *will* delete your Scratch files eventually.
:::

## Access Methods
We don't sit in front of the supercomputer. We access it remotely using one of three methods:

1.  **Open OnDemand (Browser):** The easiest way. Provides a graphical interface, file browser, and Jupyter Notebooks directly in Chrome/Firefox.
2.  **SSH (Terminal):** The "Pro" way. Fast, low bandwidth, text-only.
3.  **Virtual Desktop:** Streams a Linux desktop to your browser if you need to open GUI apps (like VMD or PyMOL).

:::{tip} Hands-on Guide available!
We have detailed step-by-step tutorials for accessing the cluster.
**Go to [Tutorial 2: Accessing Anvil via Terminal and OpenOnDemand](../tutorials/02_access_anvil.md) to get started.**
:::

:::{seealso}
Official Documentation: [Anvil Access User Guide](https://www.rcac.purdue.edu/knowledge/anvil/access)
:::