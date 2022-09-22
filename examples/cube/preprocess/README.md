# Input file generation for weak-scaling test of HemoCell

## Weak-scaling analysis

The `weakscaling.py` command-line tool can generate a directory populated with
submit scripts and simulation configuration files based on provided template
files and cluster specifications. For instance

```bash
./weakscaling.py --cluster archer2 --slurm ../templates/archer.job -v
```

populates `./batch` with a series of `weak-scaling-*` directories matching the
provided node counts set in `weakscaling.py`. The script is not fully automatic,
i.e. setting node ranges, simulation parameters, and other configurations are to
be set manually. However, after that, all template replacements and POS-file
generations are performed automatically.

Ultimately all `batch/weak-scaling-*/slurm.job` scripts can be submitted as jobs
to the compute cluster.


## Post-processing weak-scaling analysis

The `postprocessing.py` script extracts `iterate` times from timing files
generated by `HemoCell`. It handles both the traditional space-based indented
file format as well as JSON-based timing logs.