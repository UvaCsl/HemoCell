#!/usr/bin/env python3
import click
import os
import pathlib
import subprocess

from analysis import WeakScaling
from clusters import cluster_from_string, Archer2
from layout import Hexahedral
from triplet import AtomicBlock, Decomposition


def update_template_file(infile, outfile, parameters, pattern="{{{{{}}}}}"):
    """String replace patterns wrapped in three curly braces."""
    template = pathlib.Path(infile).read_text()
    for k, v in parameters.items():
        template = template.replace(pattern.format(k), str(v))
    pathlib.Path(outfile).write_text(template)


@click.command()
@click.option('--template', default="../templates/config-template.xml",
              type=click.Path(exists=True, readable=True),
              help="""Templated config.xml file.""")
@click.option('--slurm', default="../templates/hemocell-template.template-job",
              type=click.Path(exists=True, readable=True),
              help="""Templated SLURM submit script.""")
@click.option('--cluster', default='snellius',
              type=click.Choice(['snellius', 'supermuc', 'lisa', 'archer2'],
                                case_sensitive=False))
@click.option('-o', '--output', default="./batch",
              type=click.Path(dir_okay=True, writable=True),
              help="""Output directory for generate jobs files.""")
@click.option('--posfile', type=click.Path(exists=True, readable=True),
              help="""Template `RBC.pos` file used for each atomic block.""")
@click.option('-h', '--hematocrit', type=float)
@click.option('-v', '--verbose', is_flag=True)
def cli(template, slurm, cluster, output, posfile, hematocrit, verbose):
    # Parse arguments.
    cluster = cluster_from_string(cluster)
    template_file = pathlib.Path(template)
    submit_file = pathlib.Path(slurm)
    batch_directory = pathlib.Path(output)

    # Manual specification: node ranges, atomic block sizes, etc.
    nodes = (1, 4096)
    atomic_block = AtomicBlock(25, 25, 25)
    job_name = 'weak-scaling-{0:03d}'
    slurm_name = 'weak-scaling-{0:03d}.slurm'

    # Simulation specific settings to be set in the XML configuration.
    job_parameters = {
        'shear-rate': 1,
        'dx': 0.5e-6,
        'tmax': 10000,
        'tmeas': 5000,
    }

    # SLURM specific settings
    job_configuration = {
        'cpu-hour-estimate': '00:20:00',
        'batch-directory': str(batch_directory.name),
        'environment-script': cluster.environment_script,
        'energy-collection' : cluster.energy_collection,
        'ncpupernode' : cluster.cpu_per_node,
    }

    pos_file = pathlib.Path(posfile) if posfile else posfile
    domain = Decomposition(cluster.cpu_per_node)
    pattern = Hexahedral(*nodes)

    if verbose:
        click.echo(f"Hematocrit: {hematocrit}")
        click.echo(f"Posfile: {posfile}")
        click.echo(f"Config: {template_file}")
        click.echo(f"SLURM: {submit_file}")

    for i, analysis in enumerate(WeakScaling(domain, atomic_block, pattern)):
        if verbose:
            click.echo(f'{job_name.format(i)} {analysis}')

        nx, ny, nz = analysis.cells
        job_parameters.update({'nx': nx, 'ny': ny, 'nz': nz})

        casedir = batch_directory.joinpath(job_name.format(i))
        os.makedirs(casedir, exist_ok=True)
        outfile = casedir.joinpath('config.xml')

        # propagate the updates into the template files
        update_template_file(template_file, outfile, job_parameters)

        # pattern packing
        if not pos_file and not hematocrit:
            outfile = casedir.joinpath('RBC.pos')
            analysis.cell_mesh_grid(outfile)

        # cell packer packing with target hematocrit
        if not pos_file and hematocrit:
            packer_cmd = [
                "python",
                "shiftPack.py",
                "20",
                *map(str, analysis.dimensions),
                "--allowRotate",
                "--rbc",
                str(analysis.required_particles(hematocrit)),
                "--packer",
                "/home/max/work/packCells/packCells",
            ]

            if verbose:
                click.echo(" ".join(packer_cmd))

            subprocess.run(packer_cmd, check=True)
            os.rename('RBC.pos', casedir.joinpath('RBC.pos'))

        if posfile:
            # Use a template file
            outfile = casedir.joinpath('RBC.pos')
            analysis.cell_positions(posfile, outfile)

        job_configuration.update({
            'partition': cluster.partition(analysis.node_count),
            'ncpu': analysis.cpu_count,
            'nnodes' : analysis.node_count,
            'jobname': job_name.format(i),
            'slurm-file': slurm_name.format(i),
        })

        if isinstance(cluster, Archer2):
            job_configuration.update({
                'quality-of-service': cluster.qos(analysis.node_count),
                'nodes': analysis.node_count,
            })

        outfile = casedir.joinpath('slurm.job')
        update_template_file(submit_file, outfile, job_configuration)


if __name__ == "__main__":
    cli()
