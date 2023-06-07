# Welcome to MSA

MSA is short for *Matrix Structural Analysis*, and referes to the analysis of skeletal structures using the direct stiffness method.

Python package for the direct stiffness method of matrix analysis of skeletal structures. At present a subpackage `pf` for analysis of plane frames has been implemented. Code is based on the flow charts in *Weaver and Gere* and an example is taken from *Hall and Kabaila*.

## References
1. Weaver, W. and Gere, J.M., *Matrix Analysis of Framed Structures*, 2ed., CBS Publishers, New Delhi, 1986
2. Hall, A. and Kabaila, A.P., *Basic Concepts of Structural Analysis*, Pitman Publishing, London, 1977

## Install

### Clone Source from github
Clone the github repository and change into the directory. Examine the directory tree. You must have `git` installed on your machine.
```console
>git clone https://github.com/satish-annigeri/MSA.git
>cd MSA
```

### Create Virtual Environment
Create a virtual environment and install required packages.

On MS Windows do:
```console
>python -m venv venv
>venv\Scripts\activate.bat
>(venv)
```

On *nix systems do:
```console
$python3 -m venv venv
$source venv\bin
$(venv)
```

Update `pip`
```console
>(venv) python -m pip install -U pip
```

Install `pip-tools` and generate `requirements.txt` files.
```console
>(venv) pip install -U pip-tools
>(venv) pip-compile requirements.in dev-requirements.in
>(venv) pip-sync requirements.txt dev-requirements.txt
```

Execute the built-in examples in the package.
```console
>(venv) python -m msa
```

Study the `__main__.py` scripts inside the `msa` package and write your own scripts.

Finally deactivate the virtual environment.
```console
>(venv) deactivate
>
```

## Limitations
The package has the following limitations

1. Analyses only plane frame structures, with only plane frame elements.
2. Non-zero boundary conditions are currently not allowed. Only known zerodisplacement boundary condition is implemented.
3. Member loads must be applied as equivalent member end-forces. Equivalent member end-forces must be calculated manually as the negative values from the reactions produces in the corresponding fixed-end member.
4. Only member end-forces are calculates. Forces at points within the member are not calculated.