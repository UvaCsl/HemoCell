# Cube

A weak-scaling benchmark problem for `hemocell`. The problem defines a
rectangular domain subjected to a shear flow along the top and bottom planes of
the domain. Along the `x`-axis periodicity is enforced, whereas the front and
back planes of the rectangular domain have "bounce back" boundary conditions.

The problem runs for 10.000 iterations where every 1.000 a measurement step is
done, i.e. output is written and basic time per iteration is reported in the
logs. Additionally, the internal logger captures performance timings for various
components of the code throughout the full job.

## Usage

The compressed archive contains the configuration files and the initial
placement of the red blood cells (RBCs). To uncompress:

```bash
tar -xvf config.tar.xz
```

The configuration files are generated for various node counts (assuming 24 cores
per node), where their filenames match the total CPU count.

Nodes | CPUs | `nx`, `ny`, `nz` | configuration file
------|------|------------------|-------------------
1     | 24    | 150, 50, 50    | `config-24.xml`
2     | 48    | 150, 100, 50   | `config-48.xml`
4     | 96    | 150, 100, 100  | `config-96.xml`
8     | 192   | 300, 100, 100  | `config-192.xml`
16    | 384   | 300, 200, 100  | `config-384.xml`
32    | 768   | 300, 200, 200  | `config-768.xml`
64    | 1536  | 600, 200, 200  | `config-1536.xml`
128   | 3072  | 600, 400, 200  | `config-3072.xml`
256   | 6144  | 600, 400, 400  | `config-6144.xml`
512   | 12288 | 1200, 400, 400 | `config-12288.xml`
1024  | 24576 | 1200, 800, 400 | `config-24576.xml`
2048  | 49152 | 1200, 800, 800 | `config-49152.xml`

After compiling `hemocell` this example should be compiled as well and the
executable `cube` is present in this path. Then, selecting the configuration
file that matches the chosen hardware in the job definition, the problem can
be started as:

```bash
mpirun cube config/$config
```

A reference submit job (`SLURM`) is provided in the current directory
as `weakscaling.job`.

