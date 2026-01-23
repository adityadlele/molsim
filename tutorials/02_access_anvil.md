# Tutorial 2: Accessing Anvil

**Objective:** Log in to the Purdue Anvil supercomputer using the Browser (Open OnDemand) and the Terminal (SSH).

:::{note} Prerequisites
You must have your **ACCESS ID** and it must have been added to the class allocation by the instructor (see Tutorial 1).
:::

## Method 1: Open OnDemand (The Easiest Way)
**Open OnDemand (OOD)** is a web portal that allows you to use the supercomputer entirely through your web browser. No installation required.

:::{seealso}
Official Doc: [Anvil Open OnDemand Guide](https://www.rcac.purdue.edu/knowledge/anvil/access/ondemand)
:::

1.  **Navigate:** Go to [ondemand.anvil.rcac.purdue.edu](https://ondemand.anvil.rcac.purdue.edu).
2.  **Login:**
    * **Username:** Your ACCESS ID.
    * **Password:** Your ACCESS/XSEDE password.
    * **2FA:** Approve the Duo push notification.

---

## Method 2: The Virtual Desktop (For GUI Apps)
Sometimes you need to see a graphical interface (e.g., to visualize a molecule using VMD or PyMOL).

:::{seealso}
Official Doc: [Anvil Interactive Apps Guide](https://www.rcac.purdue.edu/knowledge/anvil/access/ondemand#interactive_apps)
:::

1.  On the OOD Dashboard, click **Interactive Apps** > **Anvil Desktop**.
2.  **Settings:** Set Account to `compute` and time to `1` hour.
3.  Click **Launch** and wait for the "Launch Anvil Desktop" button to appear.

---

## Method 3: SSH (The "Pro" Way)
**SSH (Secure Shell)** is how professionals access clusters. It is required for tools like VS Code Remote.

:::{seealso}
Official Doc: [Anvil SSH Access Guide](https://www.rcac.purdue.edu/knowledge/anvil/access/ssh)
:::

:::{danger} IMPORTANT: Password Restrictions
**Anvil does NOT support standard password login via terminal.**
If you try to log in via SSH without a key, your password will be rejected.
You **must** follow the "Advanced Setup" below to use SSH.
:::

### 🔑 Advanced Setup: Configuring SSH Keys
This process creates a digital "handshake" so Anvil recognizes your laptop.



#### Step A: Generate a Key (On Your Laptop)
Open your laptop terminal and run:
`ssh-keygen -t ed25519`
* Press **Enter** through all prompts.
* This creates a public key file ending in `.pub`.

#### Step B: Get your Public Key
Print the key to your screen so you can copy it.

**Mac/Linux:** `cat ~/.ssh/id_ed25519.pub`
**Windows:** `type .ssh\id_ed25519.pub`

* **Copy the output** (it usually starts with `ssh-ed25519`).

#### Step C: Install it on Anvil (The "Cheater's Method")
1.  Log in to [Open OnDemand](https://ondemand.anvil.rcac.purdue.edu).
2.  Open **Clusters** > **Anvil Shell Access**.
3.  Run these commands:
    `mkdir -p ~/.ssh`
    `nano ~/.ssh/authorized_keys`
4.  **Paste** your copied key into the editor.
5.  Press `Ctrl+O`, `Enter` to save, then `Ctrl+X` to exit.
6.  Run: `chmod 700 ~/.ssh && chmod 600 ~/.ssh/authorized_keys`

#### Step D: Test It
Back on your laptop terminal, try to log in:
`ssh username@anvil.rcac.purdue.edu`

---

## Troubleshooting
* **"Permission denied":** Your SSH key is not correctly installed in `authorized_keys`. Use Open OnDemand in the meantime.
* **"Disk Quota Exceeded":** Your Home directory is full. Move files to your Scratch space.