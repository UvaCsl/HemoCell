# Visualise cell position files

With `pos_to_vtk` cell position files used by HemoCell can be easily visualised.
This can be useful when willing to inspect cell position files before running
simulations, for instance to verify the packing results of the `packCells`
tool.

## Usage

The full usage is reported by

```bash
pos_to_vtk --help
```

For simple visualisation of a `RBC.pos`

```bash
pos_to_vtk RBC.pos
```

If you do not want to see the visualisation directly, you can provide an output
path to which the combined mesh is written to using `-o` or `--output`. The used
format is inferred based on the given path, which as to be supported by
[`meshio`](https://github.com/nschloe/meshio). For efficient file sizes and file
write speeds, using `XDMF` is suggested.

```bash
pos_to_vtk RBC.pos --output test.xdmf
```

Combining multiple `RBC.pos` files

```bash
pos_to_vtk RBC_1.pos RBC_2.pos
```

Multiple `RBC.pos` files can be combined with multiple `PLT.pos` files, where
the latter do require a separate flag per file

```bash
pos_to_vtk RBC_*.pos --plt PLT_1.pos --plt PLT_2.pos
```

## Installation

It is recommended to install the package within a [virtual
environment](https://docs.python.org/3/tutorial/venv.html). After activating
your virtual environment, run the following to install `pos_to_vtk` and its
dependencies:

```bash
pip install .

# or when you aim to edit the source code of `pos_to_vtk`
pip install --editable .
```
