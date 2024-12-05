# Python API

The FiniteFlow program comes with a Python 3 API, built on top of its C API.  This can also be used as a template for interfacing FiniteFlow to other languages or frameworks.

At the moment it is still a work in progress and the API is not stable.


## Installation

The Python API uses the [CFFI](https://cffi.readthedocs.io/en/stable/) package, which needs to be installed first.  Detailed instructions can be found [here](https://cffi.readthedocs.io/en/stable/installation.html).  On Linux systems, the following command should be enough in most cases
```
pip3 install cffi --user
```
where the option `--user` can be omitted for a global installation, if you have root privileges.  On MacOS `libffi` needs to be installed first, as explained [here](https://cffi.readthedocs.io/en/stable/installation.html#macos-x).

Moreover, FiniteFlow should also be installed before moving to the next step.

The API is then compiled using
```
cd /path/to/finiteflow/pythonapi
python3 fflow_build.py
```

After this, consider uptading your `$PYTHONPATH` by adding the following to your `~/.bashrc`, `~/.zshrc` or `~/.bash_profile`
```
export PYTHONPATH=$PYTHONPATH:/path/to/finiteflow/pythonapi
```


## Usage
The installation process described above, creates and installs the `fflow` package, which can be imported and used as any Python program.

The file `tests.py` contains several examples of usage of FiniteFlow from Python 3.

More examples will be added, once the API stabilizes.
