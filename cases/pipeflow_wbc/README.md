# pipeflow_wbc

Pipe flow with RBCs, PLTs and one WBC, for testing the numerical stability of
the cell models at high wall shear rates.

* Pipe: 20 um diameter, periodic over 61.5 um (`tube.stl`, the pipeflow tube
  stretched 1.5x along its axis; 40 nodes across at dx = 0.5 um).
* Cells: 40 RBCs (tube haematocrit about 19 %), 6 PLTs and one 8 um WBC on
  the axis at mid-length (`WBC.pos`). RBC and PLT materials are those of
  `examples/pipeflow`; `WBC.xml` is described in its comment.
* Driving: `<wallShearRate>` sets the pressure gradient G of cell-free
  Poiseuille flow with that wall shear rate, G = 2 mu gamma_w / R. In steady
  state the wall shear stress G R / 2 follows from the force balance on the
  pipe, cells or not, so with a cell-free layer at the wall the wall shear
  rate stays close to the nominal one.
* Run length: `<wallStrain>` (wall shear rate x time), at least
  `<minIterations>`, at most `<tmax>` iterations.

## Stability checks

Every `<tcheck>` iterations (a multiple of `<stepMaterialEvery>`, so the
vertex forces are those of the last material update, before the
`FORCE_LIMIT` cap) the driver reassembles all cells and writes a row to
`<outputDirectory>/stability.csv`:

* fluid: maximum lattice velocity, mean velocity [m/s];
* per cell type: number of complete cells, largest vertex force [pN],
  fraction of vertices above `FORCE_LIMIT`, largest and smallest membrane
  link strain (L/L0 - 1), largest relative triangle area and cell volume
  change, largest elongation and compression (vertex distance from the
  centroid relative to the reference mesh);
* WBC: distance of its centroid from the axis [um] and mean x velocity [m/s].

The run stops with exit code 2 and `RESULT UNSTABLE` in the log when a
position, force or velocity is not finite, a cell is lost (HemoCell deletes a
cell when a vertex enters the wall), the lattice velocity exceeds 0.3, or a
cell passes the singular points of its area (30 %) or volume (10 %) terms, or
a link strain of 2 (limits in the optional `<stability>` section). Otherwise
it ends with `RESULT STABLE` and the extreme values over the run.

## Running

```bash
cmake --build build --target pipeflow_wbc --parallel 4
cd cases/pipeflow_wbc
mpirun -n 4 ./pipeflow_wbc config.xml
```

## Initial positions

The RBC and PLT positions were placed by random sequential addition in the
pipe: random position and orientation, accepted when every vertex of the
cell's mesh (rotated as by the HemoCell loader) is at least 1 um inside the
wall, inside the periodic x range, clear of the WBC, and outside the
bounding ellipsoids of the cells placed before. Box packings from
`tools/packCells` lose most cells at the curved wall of so narrow a pipe.
