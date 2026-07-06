# Installation
<!-- To use the software, just type
```pip install pyTEM` -->

## Staying updated

To use the package and keep updated, just

1. `git clone` the repository and go to its location with bash, PowerShell etc.
2. install the code editable by `pip install -e .`
3. Update at any later time by `git pull` (no further step needed)

If you use UV, the simplest way is to create a virtual environment by typing `uv run` in the project folder which will care for all
dependencies (like empymod and pyGIMLi).
As a result, a virtual environment `.venv` is created
in the main folder and if you open the folder in VSCode it is chosen as default environment.
For use with VSCode, you should also `pip install ipykernel`.
