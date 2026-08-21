# Installation
<!-- To use the software, just type
```pip install pyTEM` -->

## Staying updated

The package will be installed like a simple `pip install pyTEMx`, only that
the name `pytem` is already reserved at PyPI.org.

The package is hosted at <https://github.com/TUBAF-EM/pyTEM>

To use the package and keep updated:

1. Clone the source code
`git clone https://github.com/TUBAF-EM/pyTEM.git`
2. go to its location with bash, PowerShell etc.
3. install the code editable by `pip install -e .`
4. Update at any later time by `git pull`

If you use UV, the simplest way is to create a virtual environment by typing
`uv run` in the project folder which will care for all dependencies (like
empymod or pyGIMLi).
As a result, a virtual environment `.venv` is created
in the main folder and if you open the folder in VSCode it is chosen as default environment.
For use with VSCode, you should also `pip install ipykernel`.
