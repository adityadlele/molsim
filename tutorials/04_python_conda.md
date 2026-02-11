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

If you want to create a *custom* environment for your own separate research projects, follow these steps:

1. **Create the environment**:
```bash
conda create -n my_md_env python=3.10 numpy matplotlib

```


2. **Activate it**:
```bash
conda activate my_md_env

```


3. **Install additional packages**:
```bash
conda install -c conda-forge mdanalysis

```


4. **Connect it to Jupyter**:
If you want your custom environment to show up in the Jupyter dropdown menu, you must install the `ipykernel` package:
```bash
conda install ipykernel
python -m ipykernel install --user --name my_md_env --display-name "Python (My_Custom_Env)"

```



---

## 📚 Resources

* [Official Conda Cheat Sheet](https://docs.conda.io/projects/conda/en/latest/user-guide/cheatsheet.html)
* [Anvil User Guide: Python](https://www.rcac.purdue.edu/knowledge/anvil/run/python)
