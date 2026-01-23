# Tutorial 3: Linux Survival Guide

**Objective:** Master the command-line skills required to navigate Anvil, manage simulation data, and troubleshoot runs.

In HPC, there is no "File Explorer." Everything is text-based. Mastering these commands is the difference between a 5-minute task and a 5-hour struggle.

---

## 1. Navigating the File System


* `pwd` (Print Working Directory): Shows your current "address."
* `ls` (List): Shows files in the current folder.
  * `ls -l`: Shows detailed info (permissions, size, date).
  * `ls -h`: Makes file sizes readable (e.g., 2G instead of 2147483648).
* `cd <folder>` (Change Directory): Move into a folder.
  * `cd ..` : Go up one level.
  * `cd ~` : Jump to Home.
  * `cd $SCRATCH` : Jump to your high-speed work area.

---

## 2. The Power of Wildcards (*)
You will often have 100s of files (e.g., `frame_1.pdb`, `frame_2.pdb`). You don't want to move them one by one.
* `ls *.pdb` : Lists **only** files ending in .pdb.
* `rm frame_*.log` : Deletes every log file starting with "frame_".
* `cp *.in ../backup/` : Copies all input files to another folder.

---

## 3. Searching and Inspecting Files
Molecular dynamics logs can be millions of lines long. Do not use `cat` on them!

* `grep "ENERGY" production.log` : Searches for the word "ENERGY" in a file and prints those lines.
* `less production.log` : Opens the file for scrolling. 
  * Use **Arrow Keys** to move.
  * Press **/** to search for a word.
  * Press **q** to quit.
* `tail -f production.log` : **"Follow"** a file. This shows you the simulation progress in real-time as the software writes to the log.

---

## 4. File Permissions (chmod)
If you see an "Access Denied" error, check your permissions.
* `chmod +x script.sh` : Makes a script "executable" so you can run it.
* `chmod 700 private_folder` : Makes a folder invisible to everyone but you.

---

## 5. Anvil Specifics: Modules


Software on Anvil is "hidden" until you ask for it. This prevents version conflicts.
* `module list` : See what is currently active.
* `module spider gromacs` : Search for GROMACS versions.
* `module load gromacs/2023.2` : Load a specific version.
* `module purge` : Clear everything (useful if your environment gets messy).

---

## 6. Pro-Tips for Speed
1. **Tab Completion:** Type `cd Sc` and hit **Tab**. It will auto-fill `cd Scratch/`. This prevents typos!
2. **Command History:** Use the **Up/Down arrows** to scroll through previous commands.
3. **The Manual:** Not sure what a command does? Type `man ls` to see the manual.

---

## 📚 External Learning Resources
If you want to dive deeper into Linux, these are the best free resources available:

* **[Linux Journey](https://linuxjourney.com/)**: Probably the best site for beginners. It’s organized by level (Grasshopper to Wizard).
* **[The Missing Semester of Your CS Education](https://missing.csail.mit.edu/2020/shell/)**: A fantastic course by MIT that covers the "shell" tools that most classes skip.
* **[Software Carpentry: The Unix Shell](https://swcarpentry.github.io/shell-novice/)**: A classic guide designed specifically for researchers and scientists.
* **[Ryans Tutorials](https://ryanstutorials.net/linuxtutorial/)**: Very easy-to-follow practical examples.

---

## Practice Exercise
1. Go to your scratch space: `cd $SCRATCH`
2. Create a "Lab1" folder: `mkdir Lab1`
3. Enter it: `cd Lab1`
4. Create three dummy files: `touch test1.txt test2.txt log.pdb`
5. List only the .txt files: `ls *.txt`
6. Check your storage usage: `myquota`