"# DSE-VS" 
This is the repository for the DSE course,  of group 10.

#### **VS Code Extensions**
- autoDocstring
- Jupyter
- GitLens
- Rainbow CSV
- vscode-pdf

####  **Setting Up the Virtual Environment**
1. Ensure that you have `conda` installed and open anaconda prompt.

2. First, we will change the conda's package-dependency solver to `libmamba`, as it is much faster. Let's first install it by activating the base environment and running

```bash

conda install conda-libmamba-solver

```

3. Deactivate the base environment and change the solver:

```bash

conda config --set solver libmamba

```

4. Now its time to install the actual libraries. Navigate the the project directory (with `cd`) and write:

```bash

conda env create --file environment.yml

```

5. Next, to ensure that python 'sees' the .py files in subfolders within project directory, activate the `(FILE DIRECTORY)` environment and then run:

```bash

conda develop .

```

6. (Optional) Install **poliastro** (good for orbital mechanics) in the X environment by activating it and running:

```bash

conda activate SVV-FD

conda install poliastro 
  
  

```

7. (Optional) If we add more packages, then run again:

```bash

conda env update --file environment.yml

```

<!-- A solver in the context of package management (like Conda) is a tool or algorithm responsible for resolving dependencies between software packages. Its primary job is to determine which versions of packages and their dependencies can coexist in an environment without conflicts. -->
