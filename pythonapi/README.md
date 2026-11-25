# Python API

The FiniteFlow program comes with a Python 3 API, built on top of its C API.  This can also be used as a template for interfacing FiniteFlow to other languages or frameworks.

At the moment, it is still a work in progress and the API is not stable.


## Installation from source

### Dependencies

The Python API uses the [CFFI](https://cffi.readthedocs.io/en/stable/) package, which needs to be installed first.  Detailed instructions can be found [here](https://cffi.readthedocs.io/en/stable/installation.html).  On Linux systems, the following command should be enough in most cases
```
pip3 install cffi --user
```
where the option `--user` can be omitted for a global installation, if you have root privileges.  On MacOS[^1] `libffi` needs to be installed first, as explained [here](https://cffi.readthedocs.io/en/stable/installation.html#macos-x).


### The Python API

Once the dependencies are installed, you can compile the Python 3 API as part of the regular installation procedure of FiniteFlow, simply by adding the option `-DFFLOW_PYTHON=1` to the `cmake` command.

Python modules will be installed in a default directory, which is suited for most use cases (it does *not* require admin privileges).  You may still optionally customize the installation directory of Python modules, using `-DFFLOW_PYTHON_INSTALL=/path/to/installation/dir`.  If you do so, consider updating your `$PYTHONPATH` by adding the following to your `~/.bashrc`, `~/.zshrc` or `~/.bash_profile`
```
export PYTHONPATH=$PYTHONPATH:/path/to/installation/dir
```


## Installation from `pip`

If you only need the Python interface of FiniteFlow, you can simply install it with the command `pip install finiteflow --user`.  See the instruction in the [main README file](../README.md) for more details.


## Usage

The installation process described above, creates and installs the `fflow` package, which can be imported and used as any Python program.

Tutorials using this interface can be found in the [tutorials/](tutorials/README.md) subdirectory.  Additional applications can be found in [`tests.py`](tests.py) and [`test_solver.py`](test_solver.py).


[^1]: When using the `python3` executable provided by Homebrew, CFFI can be installed with `brew install python-setuptools && brew install cffi`
