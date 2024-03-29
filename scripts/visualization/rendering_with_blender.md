# Rendering with Blender

Advanced visualisations of HemoCell output can be created using rendering, for
instance using Blender [1]. To enable such visualisations, the
`convert_xmf_to_x3d.py` helper script allows to convert the typical `*.xmf`
output generated by HemoCell to extensible 3D graphics (X3D) format, supported
by various computer graphics applications, such as Blender.

Typical rendering processes are iterative processes, arranging the views,
lighting, camera positions etc. However, if desired, this process can be
automated, e.g. for when processing large number of simulation time steps. For
this, the Blender Python API can be used, as illustrated in
`render_oneCellShear.py`.

## Requirements

- `pvpython`: the Python Paraview [2] API, typically comes with an installation
  of Paraview on the system.
- Blender: for rendering the X3D geometries

## One cell shear example

First run the simulation and generate the XMF output

```bash
# from hemocell/examples/oneCellShear
mpirun -n 1 oneCellShear config.xml
../../scripts/batchPostProcess.sh
```

Then convert the XMF files to X3D format

```bash
# from hemocell/examples/oneCellShear
../../scripts/visualisation/convert_xmf_to_x3d tmp/
```

An example Blender pipeline is then evaluated as
```bash
# from hemocell/examples/oneCellShear
cd ../../scripts/visualisation && blender -b -P render_oneCellShear.py
```

[1] https://www.blender.org/
[2] https://www.paraview.org/
