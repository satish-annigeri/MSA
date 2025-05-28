
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
