# Cube

A weak-scaling benchmark problem for `hemocell`. The problem defines a
rectangular domain subjected to a shear flow along the top and bottom planes of
the domain. Along the `x`-axis periodicity is enforced, whereas the front and
back planes of the rectangular domain have "bounce back" boundary conditions.

The problem runs for a preset number of iterations with a predefined output frequency.
Basic time per iteration is reported in the logs. Additionally, the internal logger captures performance timings for various components of the code throughout the full job.

## Usage

The input files including the cell positions and the job scripts can be generated by the python script included in the `preprocess` folder. For further instruction on use please refer to the `README.md` of the preprocessor script.

