#!/usr/bin/env python3
"""Smoke run of the examples and cases: a few steps of each driver with output
and a checkpoint, on a copy of its directory.

Every built driver (an executable named after its directory under examples/ or
cases/) is copied to a temporary directory, its configuration is shortened
(tmax, tmeas, tcheckpoint, tcsv, warmup; see SHORT_RUN) and it is run with
mpirun. A run passes if it ends without an error, writes HDF5 output and a
checkpoint (checkpoint.xml) if the driver calls writeOutput() and
saveCheckPoint() respectively, and its log shows no NaN. Drivers that need
inputs not in the repository, or that are too large for a short run, are
listed in SKIP with the reason; PREPARE lists inputs put into place first.

Usage (from the repository root, after building hemocell_examples and
hemocell_cases):
    scripts/ci/smoke_drivers.py [-n RANKS] [-t TIMEOUT_S] [--keep] [--slow] [DRIVER ...]
DRIVER is a directory such as examples/pipeflow; default: all built drivers
except those in SLOW (unless --slow is given).
Exit status 1 if any run fails.
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET

# Settings that bound the run length, set in every configuration that has them.
SHORT_RUN = {'tmax': 20, 'tmeas': 10, 'tcheckpoint': 20, 'tcsv': 10, 'warmup': 0, 'tbalance': 1000000,
             'minIterations': 0}

# Configuration file per driver directory (default: config.xml). Several
# entries run the driver once per configuration.
CONFIGS = {
    'examples/capillary': ['wedge_config.xml', 'bifurcation_config.xml'],
}

# Inputs to put into the driver's directory before the run: a glob pattern
# (relative to that directory) of files to copy there, or a (file name,
# contents) pair for a file to write.
PREPARE = {
    'examples/microcontraction': ['initial_states/Ht30/*.pos'],
    'cases/stent_strut': ['cell_pos_files/reference/*.pos'],
    # No positions in the repository: a run without cells.
    'cases/stentflow': [('RBC.pos', '0\n'), ('PLT.pos', '0\n')],
}

# Drivers that cannot run in a short smoke test, with the reason.
SKIP = {
    'examples/cube': 'needs a generated configuration (examples/cube/preprocess)',
    'cases/performance_testing': 'performance configurations in configs/, not a smoke test',
}

# Drivers that take minutes even for a short run, run only with --slow.
SLOW = {
    'examples/microcontraction': '13000 RBCs',
    'cases/stent_strut': 'domain of 1000 lattice nodes along the reference direction',
    'cases/unbounded': '72701 RBCs in a 256^3 domain',
}


def drivers(root):
    found = []
    for top in ('examples', 'cases'):
        for name in sorted(os.listdir(os.path.join(root, top))):
            directory = os.path.join(top, name)
            if os.access(os.path.join(root, directory, name), os.X_OK):
                found.append(directory)
    return found


def shorten(path, output):
    """Shorten the run and set <parameters><outputDirectory> (comments are dropped).
    Return the checkpoint directory, relative to the output directory."""
    tree = ET.parse(path)
    root = tree.getroot()
    for element in root.iter():
        if element.tag in SHORT_RUN:
            element.text = ' %s ' % SHORT_RUN[element.tag]
    # The top-level <parameters> (a pre-inlet has its own).
    parameters = root.find('parameters')
    if parameters is None:
        parameters = ET.SubElement(root, 'parameters')
    directory = parameters.find('outputDirectory')
    if directory is None:
        directory = ET.SubElement(parameters, 'outputDirectory')
    directory.text = output
    tree.write(path, xml_declaration=True, encoding='utf-8')
    # Where the checkpoints go, relative to the output directory.
    checkpoint = parameters.find('checkpointDirectory')
    return 'checkpoint' if checkpoint is None else checkpoint.text.strip()


def calls(directory, function):
    """Whether the driver's sources call a function."""
    for source in glob.glob(os.path.join(directory, '*.cpp')):
        with open(source, errors='replace') as f:
            if re.search(r'\b%s\s*\(' % function, f.read()):
                return True
    return False


def run(root, directory, config, ranks, timeout, keep):
    name = os.path.basename(directory)
    work = tempfile.mkdtemp(prefix='hemocell_smoke_%s_' % name)
    # The inputs (configuration, cell files, positions, geometry) are relative
    # to the driver's directory; earlier output there is not copied.
    shutil.copytree(os.path.join(root, directory), os.path.join(work, name), symlinks=True,
                    ignore=shutil.ignore_patterns('tmp*', 'output*', 'smoke_out*'))
    cwd = os.path.join(work, name)
    for item in PREPARE.get(directory, []):
        if isinstance(item, tuple):
            with open(os.path.join(cwd, item[0]), 'w') as f:
                f.write(item[1])
        else:
            for path in glob.glob(os.path.join(cwd, item)):
                shutil.copy(path, cwd)
    checkpoint = shorten(os.path.join(cwd, config), 'smoke_out')
    command = ['mpirun', '-n', str(ranks), './' + name, config]
    start = time.time()
    try:
        result = subprocess.run(command, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                timeout=timeout)
        log, status = result.stdout.decode(errors='replace'), result.returncode
    except subprocess.TimeoutExpired as e:
        log, status = (e.stdout or b'').decode(errors='replace'), 'timeout'
    seconds = time.time() - start
    problems = []
    if status != 0:
        problems.append('exit status %s' % status)
    out = os.path.join(cwd, 'smoke_out')
    if calls(cwd, 'saveCheckPoint') and not os.path.isfile(os.path.join(out, checkpoint, 'checkpoint.xml')):
        problems.append('no checkpoint')
    hdf5 = os.path.join(out, 'hdf5')
    if calls(cwd, 'writeOutput') and not (os.path.isdir(hdf5) and any(files for _, _, files in os.walk(hdf5))):
        problems.append('no HDF5 output')
    if re.search(r'\bnan\b', log, re.IGNORECASE):
        problems.append('NaN in the log')
    if not keep:
        shutil.rmtree(work, ignore_errors=True)
    return problems, seconds, log, work


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-n', '--ranks', type=int, default=2)
    parser.add_argument('-t', '--timeout', type=int, default=600, help='per run [s]')
    parser.add_argument('--keep', action='store_true', help='keep the run directories')
    parser.add_argument('--slow', action='store_true', help='also run the drivers listed in SLOW')
    parser.add_argument('drivers', nargs='*')
    args = parser.parse_args()
    root = os.getcwd()
    selected = [d.rstrip('/') for d in args.drivers] or drivers(root)
    failed = 0
    for directory in selected:
        if directory in SKIP:
            print('SKIP %-45s %s' % (directory, SKIP[directory]), flush=True)
            continue
        if directory in SLOW and not (args.slow or args.drivers):
            print('SKIP %-45s slow (%s), run with --slow' % (directory, SLOW[directory]), flush=True)
            continue
        for config in CONFIGS.get(directory, ['config.xml']):
            label = directory if config == 'config.xml' else '%s (%s)' % (directory, config)
            problems, seconds, log, work = run(root, directory, config, args.ranks, args.timeout, args.keep)
            if problems:
                failed += 1
                print('FAIL %-45s %6.1f s  %s' % (label, seconds, ', '.join(problems)), flush=True)
                # The error messages (stderr, unbuffered) can precede the
                # last lines of the (buffered) standard output.
                lines = log.strip().splitlines()
                errors = [line for line in lines[:-15] if re.search(r'error|exception', line, re.IGNORECASE)]
                print('\n'.join('     | ' + line for line in errors[-5:] + lines[-15:]), flush=True)
            else:
                print('OK   %-45s %6.1f s' % (label, seconds), flush=True)
            if args.keep:
                print('     kept in %s' % work)
    print('%d failed' % failed)
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
