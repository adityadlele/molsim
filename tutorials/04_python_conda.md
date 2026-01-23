# Tutorial 4: Python, Conda, and Jupyter

**Objective:** Set up a reproducible Python environment on Anvil and launch a Jupyter Notebook using Open OnDemand.

---

## 1. Why use Conda?
On an HPC, you cannot use "pip install" to global directories because you don't have "Admin" (root) privileges. **Conda** allows you to create your own isolated "bubbles" (environments) where you can install specific versions of libraries (like MDAnalysis, NumPy, or PyTorch) without breaking the system software or other students' projects.



---

## 2. Using the Class Environment (Terminal)
To save time, a pre-configured environment has been created for this course. To use it, you must first tell the system where the class modules are located, then load the specific environment.

**To activate the class environment:**
1. Log in to the Anvil Terminal.
2. Run the following two commands:
   module use /anvil/projects/x-chm250117/etc/modules
   module load conda-env/molsimclass-py3.12.8

Once loaded, your terminal prompt should change to indicate the environment is active.

---

## 3. Running Jupyter Notebooks with the Class Environment
You don't run Jupyter from the terminal. Instead, use the web portal to get a dedicated compute node that automatically loads the class tools.

1. Go to: **https://ondemand.anvil.rcac.purdue.edu**
2. Click **Interactive Apps** > **Jupyter Notebook**.
3. **Settings:**
   - **Account:** compute (or your class allocation).
   - **Partition:** shared.
   - **Number of Hours:** 1 or 2.
4. **The "Conda Environment" Box:**
   - To use the class environment, paste the following path into the "Optional: Path to Conda Environment" field:
     /anvil/projects/x-chm250117/etc/conda-envs/molsimclass-py3.12.8
5. Click **Launch**.
6. Once the session is ready, click **Connect to Jupyter**.



---

## 4. Creating Your Own Environment (Optional)
If you want to create a custom environment for your own research projects:

1. **Create the environment**:
   conda create -n my_md_env python=3.10 numpy matplotlib

2. **Activate it**:
   conda activate my_md_env

3. **Install additional packages**:
   Once activated, you can add more tools:
   conda install -c conda-forge mdanalysis

---

## 5. Connecting your Custom Env to Jupyter (IPyKernel)
If you created a **custom** environment (not the class one) and it doesn't show up in the Jupyter "New" menu, run this **once** inside your activated environment in the terminal:

conda install ipykernel
python -m ipykernel install --user --name my_md_env --display-name "Python (My_Custom_Env)"

Now, when you refresh Jupyter, "Python (My_Custom_Env)" will appear as an option.

---

## 📚 Resources
* [Official Conda Cheat Sheet](https://docs.conda.io/projects/conda/en/latest/user-guide/cheatsheet.html)
* [Anvil Python Guide](https://www.rcac.purdue.edu/knowledge/anvil/run/python)