import numpy as np
import pathlib
from cell import Cell
from math import prod
from triplet import Triplet


class Analysis:
    def __init__(self, domain, atomic_block, pattern, lattice_unit=0.5):
        self.domain = domain
        self.atomic_block = atomic_block
        self.blocks = pattern
        self.lattice_unit = lattice_unit

    def __repr__(self):
        blocks = f'{self.blocks.x}x{self.blocks.y}x{self.blocks.z}'
        domain = f'{self.domain.x}x{self.domain.y}x{self.domain.z}'
        return (
            f"node: {self.node_count:>5}, "
            f"cpu: {self.cpu_count:>8}, "
            f"domain: {domain}, "
            f"blocks: {blocks:>8}, "
            f"LBM: {str(self.cells):>14}, "
            f"dimensions: {str(self.dimensions):>20}"
        )

    @property
    def cells(self):
        """Number of lattice Boltzmann cells in the domain."""
        return self.domain * self.atomic_block * self.blocks

    @property
    def dimensions(self):
        """The physical dimensions of the domain in micron."""
        return self.cells * self.lattice_unit

    @property
    def cell_count(self):
        """The total count of lattice Boltzmann cells in the domain."""
        return prod(self.cells)

    @property
    def node_count(self):
        """The number of cluster nodes used."""
        return prod(self.blocks)

    @property
    def cpu_count(self):
        """The total number of CPUs used for the domain."""
        return prod(self.domain) * self.node_count

    def offsets(self):
        """The required offsets to shift a templated cell position file."""
        atomic_block_size = self.domain * self.atomic_block * self.lattice_unit
        for z in range(self.blocks.z):
            for y in range(self.blocks.y):
                for x in range(self.blocks.x):
                    yield atomic_block_size * Triplet(x, y, z)

    def cell_mesh_grid(self, outfile):
        """Generates a cell packing along a fixed mesh-grid orientation.

        This creates a cell packing from a given cell count in x, y, z
        direction for the smallest atomic block. The pattern of those cells is
        then repeated throughout the domain. This has the advantage of
        providing a completely constant hematocrit values across domain sizes,
        with the disadvantage of little cells across node boundaries with
        respect to packing obtained from the cell packer.
        """
        cell_diameters = Triplet(8.4, 4.4, 8.4)
        cell_angles = Triplet(0, 0, 0)

        # The pattern to be fitted within the domain
        pattern = Triplet(5, 9, 5) * self.blocks

        def pad(ncells, domain_length, cell_diameter):
            """Required padding to fill the domain.

            Ensures the specified pattern is centered within the domain by
            introducing a slight packing at the start/end.
            """
            return domain_length/ncells/2 - cell_diameter/2

        padding = Triplet(*map(pad, pattern, self.dimensions, cell_diameters))

        # Assert the padded domain remains within the specified dimensions
        padded_size = pattern * (padding * 2 + cell_diameters)
        assert all((a <= b for (a, b) in zip(padded_size, self.dimensions)))

        # The resulting `x, y, z` cell positions
        xr, yr, zr = [((p+D/2), d, (2*p+D)) for (p, d, D)
                      in zip(padding, self.dimensions, cell_diameters)]

        # The final cell positions are found by the product from the padded
        # ranges that were obtained previously.
        cells = []
        for z in np.arange(*zr):
            for y in np.arange(*yr):
                for x in np.arange(*xr):
                    cells.append(Cell(x, y, z, *cell_angles))

        with open(outfile, 'w') as file:
            file.write(f'{len(cells)}\n')
            for cell in cells:
                file.write(f'{cell}\n')

    def cell_positions(self, pos_file, outfile):
        """Use a cell position file as template to pack each node's domain.

        The given POS-file is applied as template for all atomic blocks and
        repeated across the full simulation domain. This allows to generate
        random packings for the smallest domain and repeat these consistently,
        however, has the downside of a lack of cells across node boundaries.
        """
        pos_file = pathlib.Path(pos_file).read_text()
        cells = (Cell.from_str(s) for s in pos_file.splitlines())
        cells = list(filter(None, cells))

        cell_count = prod(self.blocks) * len(cells)
        with open(outfile, 'w') as file:
            file.write(f'{cell_count}\n')
            for cell in cells:
                for offset in self.offsets():
                    file.write(f'{cell.translate(*offset)}\n')

    def required_particles(self, hematocrit):
        """Estimated particle count to achieve desired hematocrit."""
        volume = prod(self.cells * self.lattice_unit)
        particle_count = (hematocrit * volume)/90
        return int(round(particle_count))


class WeakScaling:
    """Yields analysis instances for the full weak-scaling range."""

    def __init__(self, decomposition, atomic_block, layout):
        self.dec = decomposition
        self.atomic_block = atomic_block
        self.layout = layout

    def __iter__(self):
        for layout in self.layout:
            yield Analysis(self.dec, self.atomic_block, layout)
