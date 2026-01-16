# Moltemplate: four-layer Ag(111) slab workflow

This repository already uses ASE to build metal slabs (see `Pt111.py`/`pt3_2.py`).
If you want to build a **four-layer Ag(111) surface** for use in moltemplate, an
easy and reliable workflow is:

1. **Generate the Ag(111) slab with ASE** (4 layers, desired in-plane size, and
   optional vacuum).
2. **Export to a LAMMPS data file**.
3. **Import that data file from moltemplate** and compose your full system in
   `system.lt`.

Below is a minimal example that matches that flow and can be adjusted for your
cell size, lattice constant, and vacuum.

---

## 1) Generate a 4-layer Ag(111) slab with ASE

Create a Python script (example: `Ag111_4layer.py`) and run it to export a
LAMMPS data file.

```python
from ase.build import fcc111
from ase.io import write

# FCC Ag lattice constant (angstrom). Adjust if you have a preferred value.
a = 4.09

# size=(nx, ny, nlayers) -> here nlayers=4
# vacuum adds empty space above the slab along z
ag111 = fcc111('Ag', size=(3, 3, 4), a=a, vacuum=12.0, orthogonal=True)

# Write LAMMPS data file (atomic style)
write('ag111_4layer.data', ag111, format='lammps-data', atom_style='atomic')
```

Notes:
- Change `size=(3, 3, 4)` to the in-plane replication you want.
- `orthogonal=True` yields an orthogonal cell (easier for moltemplate/LAMMPS).
- Increase/decrease `vacuum` depending on your simulation setup.

---

## 2) Import the slab into moltemplate

In moltemplate, you can import the LAMMPS data file and then build the rest of
your system. A minimal `system.lt` might look like this:

```lt
# system.lt

import "ag111_4layer.data"

# You can add other molecules or regions here and position them
# relative to the Ag(111) slab.
```

Then run moltemplate normally:

```bash
moltemplate.sh system.lt
```

This produces `system.data` (and optionally `system.in*`) which you can run in
LAMMPS.

---

## 3) Tips for customization

- **Constrain bottom layers**: after import, you can add `fix` commands in your
  LAMMPS input to freeze the bottom 1–2 layers.
- **Surface orientation**: `fcc111` is the standard (111) surface for fcc
  lattices; using `orthogonal=True` avoids the 60° tilt in the simulation cell.
- **Lattice constant**: update `a` if you are using a specific potential with a
  different equilibrium lattice parameter.

---

If you'd rather build the slab directly in moltemplate (without ASE), you can do
it by enumerating the Ag lattice positions manually in an `.lt` file, but the
ASE → LAMMPS data → moltemplate import workflow is usually much faster and less
error-prone for crystal slabs.
