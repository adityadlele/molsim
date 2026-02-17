# Tutorial 4: Python, Conda, and Jupyter

**Objective:** Set up the reproducible class Python environment on Anvil and launch a Jupyter Notebook.

---

## 1. Why use Conda?
On an HPC (High Performance Computer), you cannot simply `pip install` packages globally because you do not have "Admin" (root) privileges. If everyone installed different versions of NumPy or PyTorch into the system folders, everything would break.

**Conda** solves this by creating isolated "bubbles" (environments).
* **Isolation:** You can have one bubble for this class (with `numpy 1.26`) and another for your research (with `numpy 2.0`) without conflict.
* **Reproducibility:** Everyone in the class uses the exact same bubble, so code that runs for me will run for you.



---

## 2. Setting Up the Class Environment
To ensure compatibility, we have created a shared central environment on Anvil. You do not need to install anything yourself. You simply need to tell Anvil where to find our class "bubble."

### Step 1: Open a Terminal
1.  Log in to [Anvil OnDemand](https://ondemand.anvil.rcac.purdue.edu).
2.  Click **Clusters** > **Anvil Shell Access** (or open a Terminal within an existing Jupyter session).

### Step 2: Run the Setup Commands
Copy and paste the following commands into your terminal. You only need to run these **once** to set up your account.

```bash
# 1. Tell Anvil where to look for our class modules
module use /anvil/projects/x-chm250117/etc/modules

# 2. Load the specific class environment module
module load conda-env/molsimclass-py3.12.8

# 3. Register the environment as a Jupyter Kernel
conda-env-mod kernel -p /anvil/projects/x-chm250117/apps/molsimclass

```

**What did I just do?**

* `module use`: You added our class's private software folder to your path.
* `module load`: You actually loaded the software environment.
* `conda-env-mod kernel`: You created a shortcut (a "Kernel") so the Jupyter web interface knows this Python environment exists.

---

## 3. Launching Jupyter with the Class Kernel

Now that you have registered the environment, you can use it in your notebooks.

1. Go to **[Anvil OnDemand](https://ondemand.anvil.rcac.purdue.edu)**.
2. Click **Interactive Apps** > **Jupyter Notebook**.
3. **Settings:**
* **Account:** `compute` (or your specific allocation).
* **Partition:** `shared`.
* **Number of Hours:** `1` or `2`.


4. Click **Launch**.
5. Once the session is ready, click **Connect to Jupyter**.

### Selecting the Kernel

When your notebook opens, it might default to the standard "Python 3". To use the class tools:

1. Open your notebook.
2. In the top right corner, click on the kernel name (e.g., "Python 3").
3. In the dropdown menu, select: **`molsimclass`** (or `Python [conda env: molsimclass]`).

**Verification:**
Run this Python code in a cell. It should point to the project folder, not your home folder:

```python
import sys
print(sys.executable)
# Expected output: /anvil/projects/x-chm250117/apps/molsimclass/...

```

---

## 4. Creating Your Own Environment (Optional)

If you want to create a *custom* environment for your own separate research projects, we recommend using the system's `conda-env-mod` tool. Unlike standard conda commands, this tool automatically generates a **module file** (so you can easily load it later) and a **Jupyter kernel** (so it appears in your notebook dropdowns).

**Step 1: Load the Conda module**
You must load the system's base conda module to access the creation script.

```bash
module load conda

```

**Step 2: Create the environment**
Use the `conda-env-mod` script to create your environment. We recommend adding the `--jupyter` flag immediately so it sets up the kernel for you.

```bash
# Create an environment named 'my_md_env' with Jupyter support
conda-env-mod create -n my_md_env --jupyter

```

*Follow the on-screen prompts. You may need to type `y` to confirm the creation.*

**Step 3: Load your new environment**
Once created, the script will print specific instructions on how to load your new environment. You must load this module to use it. It generally follows this pattern:

```bash
# 1. Tell the system where to look for your custom modules
module use $HOME/privatemodules

# 2. Load your specific environment module
# (Note: The actual version suffix 'py3.X.X' will be displayed in the output of Step 2)
module load conda-env/my_md_env-py3.8.8

```

**Step 4: Install Additional Packages**
Once the module is loaded (Step 3), you are "inside" your environment. You can now use standard `conda` or `pip` commands to install any software you need.

* **Install via Conda (Recommended for scientific libraries):**
```bash
# Install standard libraries
conda install numpy matplotlib

# Install MDAnalysis from the conda-forge channel
conda install -c conda-forge mdanalysis

```


* **Install via Pip (If not available on Conda):**
```bash
pip install package_name

```



**Note on Jupyter:** Because you used the `--jupyter` flag in Step 2, you do **not** need to manually install `ipykernel`. Your environment will automatically appear as **"Python (My my_md_env Kernel)"** in the Jupyter interface.


**Detailed Documentation:**
For more advanced usage, troubleshooting, or cluster-specific policies, please refer to the official [Anvil User Guide: Installing Python Packages](https://www.rcac.purdue.edu/knowledge/anvil/run/examples/apps/python/packages).





## 📚 Resources

* [Official Conda Cheat Sheet](https://docs.conda.io/projects/conda/en/latest/user-guide/cheatsheet.html)
* [Anvil User Guide: Python](https://www.rcac.purdue.edu/knowledge/anvil/run/examples/apps/python/packages)
```