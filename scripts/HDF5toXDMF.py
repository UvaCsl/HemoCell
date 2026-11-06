#!/usr/bin/env python3
"""Write XDMF (version 3) descriptions of HemoCell HDF5 output.

HemoCell writes these files itself while it runs. Use this script for output of
older HemoCell versions (or to rewrite the files). It writes, next to the hdf5/
directory of an output directory, for every output name (cell type, Fluid,
CEPAC, with _PRE for a pre-inlet):

  <name>.<iteration>.xmf  one output step: a spatial collection of the blocks
  <name>.xmf              all steps: a temporal collection including the step files

Open <name>.xmf in ParaView (any of its XDMF readers) to get all steps with
their time.

Fluid blocks are written as the cells of a regular grid, one cell per lattice
node. Output of older versions holds one extra node layer around each block;
it is shown as well (neighbouring blocks overlap by those layers), and cells
were written by every block holding part of them, so they may appear several
times.

Requires h5py. Usage:

  HDF5toXDMF.py [--force] [directory ...]

A directory is an output directory (with an hdf5/ subdirectory) or a directory
whose subdirectories are output directories (e.g. a case directory with tmp/,
tmp_0/); the default is the current directory. Existing step files are kept
unless --force is given; the <name>.xmf files are always rewritten.
"""
import argparse
import os
import re
import sys
from collections import defaultdict
from xml.sax.saxutils import escape, quoteattr

import h5py

BLOCK_FILE = re.compile(r'^(?P<name>.+)\.(?P<iteration>\d+)\.p(?:gz)?\.(?P<block>\d+)\.h5$')
STEP_FILE = r'^{}\.\d+\.xmf$'
ATTRIBUTE_TYPE = {1: 'Scalar', 3: 'Vector', 6: 'Tensor6', 9: 'Tensor'}
TOPOLOGY_DATASETS = ('Triangles', 'InnerLinks')


def number(value):
    return '%.15g' % value


def hdf_item(dimensions, h5dataset, path, dataset):
    # NumberType and Precision from the stored type (e.g. "Cell Id" is a
    # 64-bit integer since HemoCell 3.0, older files store float).
    kind = h5dataset.dtype.kind
    number_type = 'Float' if kind == 'f' else ('UInt' if kind == 'u' else 'Int')
    return ('<DataItem Dimensions="%s" NumberType="%s" Precision="%d" Format="HDF">%s:/%s</DataItem>'
            % (' '.join(str(d) for d in dimensions), number_type, h5dataset.dtype.itemsize, escape(path),
               escape(dataset)))


def attribute(name, components, center, item):
    return ('        <Attribute Name=%s AttributeType="%s" Center="%s">\n          %s\n        </Attribute>\n'
            % (quoteattr(name), ATTRIBUTE_TYPE.get(components, 'Matrix'), center, item))


def structured_block(h5, path, block):
    nodes = [int(n) for n in h5.attrs['subdomainSize']]  # z, y, x
    spacing = [float(v) for v in h5.attrs['dxdydz']]
    # The attributes are single precision. The spacing is dx (SI) or 1
    # (lattice units), and cell corners lie half way between nodes, so both
    # can be restored exactly.
    dx = float(h5.attrs['dx'][0]) if 'dx' in h5.attrs else None
    spacing = [dx if dx and abs(s - dx) <= 1e-6 * dx else s for s in spacing]
    origin = [round(2 * float(o) / s) / 2 * s for o, s in zip(h5.attrs['relativePosition'], spacing)]
    out = '      <Grid Name="block %d" GridType="Uniform">\n' % block
    # One cell per node: one point more than nodes in every direction.
    out += '        <Topology TopologyType="3DCoRectMesh" Dimensions="%d %d %d"/>\n' % tuple(n + 1 for n in nodes)
    out += '        <Geometry GeometryType="ORIGIN_DXDYDZ">\n'
    for values in (origin, spacing):
        out += ('          <DataItem Dimensions="3" NumberType="Float" Precision="8" Format="XML">%s</DataItem>\n'
                % ' '.join(number(v) for v in values))
    out += '        </Geometry>\n'
    for name, dataset in h5.items():
        components = dataset.shape[3] if len(dataset.shape) == 4 else 1
        out += attribute(name, components, 'Cell', hdf_item(nodes + [components], dataset, path, name))
    out += '      </Grid>\n'
    return out


def cell_block(h5, path, block):
    if 'Position' not in h5 or len(h5['Position'].shape) != 2 or h5['Position'].shape[0] == 0:
        return ''
    vertices = h5['Position'].shape[0]
    position = hdf_item([vertices, 3], h5['Position'], path, 'Position')
    out = '      <Grid Name="block %d" GridType="Uniform">\n' % block
    triangles = h5['Triangles'].shape[0] if 'Triangles' in h5 else 0
    if triangles > 0:
        out += '        <Topology TopologyType="Triangle" NumberOfElements="%d">\n' % triangles
        out += '          %s\n        </Topology>\n' % hdf_item([triangles, 3], h5['Triangles'], path, 'Triangles')
    else:
        out += ('        <Topology TopologyType="Polyvertex" NumberOfElements="%d" NodesPerElement="1"/>\n'
                % vertices)
    out += '        <Geometry GeometryType="XYZ">\n          %s\n        </Geometry>\n' % position
    for name, dataset in h5.items():
        if name in TOPOLOGY_DATASETS or dataset.shape[0] != vertices:
            continue
        components = dataset.shape[1] if len(dataset.shape) == 2 else 1
        out += attribute(name, components, 'Node', hdf_item([vertices, components], dataset, path, name))
    out += '      </Grid>\n'
    links = h5['InnerLinks'].shape[0] if 'InnerLinks' in h5 else 0
    if links > 0:
        out += '      <Grid Name="block %d inner links" GridType="Uniform">\n' % block
        out += ('        <Topology TopologyType="Polyline" NumberOfElements="%d" NodesPerElement="2">\n'
                '          %s\n        </Topology>\n' % (links, hdf_item([links, 2], h5['InnerLinks'], path, 'InnerLinks')))
        out += '        <Geometry GeometryType="XYZ">\n          %s\n        </Geometry>\n' % position
        out += '      </Grid>\n'
    return out


def write_file(path, content):
    temporary = path + '.tmp'
    with open(temporary, 'w') as f:
        f.write(content)
    os.replace(temporary, path)


def write_step(directory, name, iteration, blocks):
    """Write <name>.<iteration>.xmf for the block files [(block id, relative path)]."""
    grids = ''
    time = None
    for block, path in sorted(blocks):
        with h5py.File(os.path.join(directory, path), 'r') as h5:
            if time is None:
                it = h5.attrs.get('iteration')
                dt = h5.attrs.get('dt')
                time = float(it[0]) * float(dt[0]) if it is not None and dt is not None else float(iteration)
            grids += structured_block(h5, path, block) if 'subdomainSize' in h5.attrs else cell_block(h5, path, block)
    content = ('<?xml version="1.0" ?>\n<Xdmf Version="3.0">\n  <Domain>\n'
               '    <Grid Name=%s GridType="Collection" CollectionType="Spatial">\n'
               '      <Time Value="%s"/>\n%s    </Grid>\n  </Domain>\n</Xdmf>\n'
               % (quoteattr(name), '%.12g' % time, grids))
    write_file(os.path.join(directory, '%s.%s.xmf' % (name, iteration)), content)


def write_series(directory, name):
    """Write <name>.xmf including every step file of name in directory."""
    pattern = re.compile(STEP_FILE.format(re.escape(name)))
    steps = sorted(f for f in os.listdir(directory) if pattern.match(f))
    includes = ''.join('      <xi:include href=%s xpointer="element(/1/1/1)"/>\n' % quoteattr(s) for s in steps)
    content = ('<?xml version="1.0" ?>\n<Xdmf Version="3.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n'
               '  <Domain>\n    <Grid Name=%s GridType="Collection" CollectionType="Temporal">\n'
               '%s    </Grid>\n  </Domain>\n</Xdmf>\n' % (quoteattr(name), includes))
    write_file(os.path.join(directory, name + '.xmf'), content)
    return len(steps)


def convert(directory, force):
    hdf5 = os.path.join(directory, 'hdf5')
    steps = defaultdict(list)  # (name, iteration) -> [(block, path relative to directory)]
    for step in sorted(os.listdir(hdf5)):
        if not os.path.isdir(os.path.join(hdf5, step)):
            continue
        for f in os.listdir(os.path.join(hdf5, step)):
            m = BLOCK_FILE.match(f)
            if m:
                steps[(m.group('name'), m.group('iteration'))].append(
                    (int(m.group('block')), 'hdf5/%s/%s' % (step, f)))
    written = defaultdict(int)
    for (name, iteration), blocks in sorted(steps.items()):
        target = os.path.join(directory, '%s.%s.xmf' % (name, iteration))
        if force or not os.path.exists(target):
            write_step(directory, name, iteration, blocks)
            written[name] += 1
    for name in sorted(set(n for n, _ in steps)):
        total = write_series(directory, name)
        print('%s: %s.xmf (%d steps, %d step files written)' % (directory, name, total, written[name]))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('directories', nargs='*', default=['.'], help='output or case directories (default: .)')
    parser.add_argument('--force', action='store_true', help='rewrite existing step files')
    args = parser.parse_args()
    found = False
    for top in args.directories:
        candidates = [top] + sorted(os.path.join(top, d) for d in os.listdir(top)
                                    if os.path.isdir(os.path.join(top, d)))
        for directory in candidates:
            if os.path.isdir(os.path.join(directory, 'hdf5')):
                convert(directory, args.force)
                found = True
    if not found:
        sys.exit('No output directory (with an hdf5/ subdirectory) found in: ' + ' '.join(args.directories))


if __name__ == '__main__':
    main()
