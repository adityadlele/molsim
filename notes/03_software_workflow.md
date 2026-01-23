# HPC Software & Workflow Logic

To use a supercomputer effectively, you must understand how software is managed. Unlike a personal computer where you "Install and Forget," HPC software is modular and environment-based.

---

## 1. Why Linux?
Supercomputers run Linux because it is:
* **Efficient:** No background updates or heavy GUI slowing down the CPU.
* **Scriptable:** You can automate 1,000 simulations using a single text script.
* **Multi-user:** It is built from the ground up to allow hundreds of people to share the same machine securely.



---

## 2. The Software Stack (Modules)
In a shared environment, two researchers might need different versions of the same software (e.g., Python 3.8 vs. Python 3.12). If we "installed" both normally, they would conflict.

**The Solution: Lmod (Environment Modules)**
Software on Anvil is stored in "containers" that are invisible by default. When you run `module load`, you are telling Linux to temporarily add a specific software's "path" to your current session.

* **Key Concept:** Modules only last as long as your current terminal session. If you log out and back in, you must load them again.

---

## 3. Python & Conda: The "Environment" Philosophy
In Molecular Dynamics, we rely on many moving parts (NumPy, MDAnalysis, RDKit, etc.). 

If you update one library, it might break another. To prevent this, we use **Conda Environments**. Think of an environment as a "project-specific toolbox."
1. You create a toolbox for "Class Lab 1."
2. You install exactly what you need.
3. You "close" the toolbox when done.



---

## 4. The Lifecycle of an MD Job
When you perform research on Anvil, you typically follow this workflow:

1. **Preparation (Login Node):** You log in via SSH or Open OnDemand. You write your scripts, clean your PDB files, and check your quotas.
2. **Environment Setup:** You load your modules (`module load conda`) and activate your environment.
3. **Execution (Compute Nodes):** You do NOT run heavy simulations on the login node. You submit a "Job" to the **Batch Scheduler (SLURM)**, which finds free compute nodes for you.
4. **Analysis:** Once the simulation is done, you use Jupyter Notebooks to analyze the results and create graphs.

---

## 📚 Key Takeaways
* **Compute is remote:** You are sending instructions to a machine in a data center.
* **Storage is tiered:** Use `$SCRATCH` for work and `$HOME` for keeps.
* **Software is modular:** Always check your `module list` and `conda env list` if something isn't working.

:::{tip} Ready to practice?
Now that you understand the "Why," head over to the **Tutorials** section to start the "How."
:::