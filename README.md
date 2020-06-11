# Spice Stuff for PyMess

This module is used to obtain information such as Mercury's position around the Sun, Carrington longitude etc.

## Usage

Firstly a few environment variables need setting:

```bash
export SPICE_KERNEL_PATH=/path/to/spice/kernels
export SPICE_OUTPUT_PATH=/path/to/output/stuff
```

Also, some dependencies need installing:

```bash
pip3 install spiceypy RecarrayTools PyFileIO numpy scipy DateTimeTools --user
```

Then we can import the module in Python 3:

```python
import PlanetSpice as ps
```

where `ps` currently contains the following submodules:

[`Sun`](PlanetSpice/Sun/README.md)

[`Mercury`](PlanetSpice/Mercury/README.md)

[`Venus`](PlanetSpice/Venus/README.md)

[`Earth`](PlanetSpice/Earth/README.md)

[`Mars`](PlanetSpice/Mercury/README.md)

 