#!/usr/bin/env python

import click
import json
import numpy as np
import pathlib


def parse(filepath):
    try:
        data = json.loads(filepath.read_bytes())
        data = [v['HemoCell']['iterate']['Total'] for v in data.values()]
    except json.JSONDecodeError:
        contents = filepath.read_text().splitlines()
        contents = filter(lambda l: 'iterate' in l, contents)
        data = list(map(lambda l: float(l.strip().split()[-1]), contents))
    return data


def efficiency(times, start=0):
    reference = times[start]
    return [reference/t for t in times]


def plot(node_count, eta, cpu_per_node, show=False, output=None):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise click.UsageError('Cannot import Matplotlib.')

    plt.plot(node_count, eta, marker='*')
    plt.ylim(0, 1.1)
    plt.xlabel(f"Node (x{cpu_per_node} CPU)")
    plt.ylabel("Efficiency")

    if output:
        plt.savefig(output)

    if show:
        plt.show()


@click.command()
@click.argument('samples', nargs=-1, type=click.Path(exists=True))
@click.argument('cpu-per-node', type=int)
@click.option('--show', 'show', is_flag=True,
              help="""Show interactive matplotlib efficiency figure.""")
@click.option('-o', '--output', 'output', type=click.Path(writable=True),
              help="""Write matpltolib figure to disk""")
@click.option('--sort', is_flag=True, default=False,
              help="""Sort list of SAMPLE file alfabetically. Take care for
              files with/without leading zeros before node counts.""")
def cli(cpu_per_node, samples, show, output, sort):
    """Post-process HemoCell timing samples.

    The TIMINGS files can be provided in the (legacy) indented format or the
    recent JSON format. This extracts iterate times only to estimate efficiency
    and performance metrics. The CPU-PER-NODE count is required to determine
    node counts for each job.
    """
    samples = list((pathlib.Path(sample).resolve() for sample in samples))
    samples = sorted(samples) if sort else samples

    # Parse data
    times = list(map(parse, samples))
    node_count = [len(d)//cpu_per_node for d in times]

    # Extract performance
    mean = list(map(np.mean, times))
    std = list(map(np.std, times))
    eta = [mean[1]/m for m in mean]
    speedup = [n*2**(i) for i, n in enumerate(eta)]

    # Dump tab-separated performance summary. Can be transformed to other
    # formats using `cut` and `tr` tools, e.g. extract only node count and
    # efficiency in comma separated list:
    #
    # ./postprocess.py ... | cut -f 2,5 | tr '\t' ','
    for i in range(len(samples)):
        print("{}\t{}\t{:.2f}\t({:.2f})\t{:.2f}\t{:.2f}".format(
            samples[i].name,
            node_count[i],
            mean[i],
            std[i],
            eta[i],
            speedup[i]
        ))

    if show or output:
        plot(node_count, eta, cpu_per_node, show, output)


if __name__ == "__main__":
    cli()
