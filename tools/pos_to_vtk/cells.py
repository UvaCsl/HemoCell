from cell import Cell
import click
from tqdm import tqdm


class Cells:
    """A collection of `cell.Cell`s.

    The collection provides some utility functions to clip the cells with
    meshes or bounding boxes and provide an hematocrit estimate for those
    volumes as well. This class is typically initialised through it `from_pos`
    static method that initialises all cells for a given position file.
    """

    def __init__(self, cells):
        self.cells = cells

    def __repr__(self):
        return f'Cells with count: {len(self)}'

    def __iter__(self):
        for cell in self.cells:
            yield cell

    def __len__(self):
        return len(self.cells)

    @staticmethod
    def from_pos(path):
        """Initialises all cells in the pos file that successfully parse.

        Any cells that fail to parse are ignored.
        """
        # Note: need to use `click.open_file` here to allow passing a dash
        # (`-`) as position file argument to indicate cell positions are given
        # through `stdin`. Using plain `open` will fail here otherwise.
        with click.open_file(path, 'r') as lines:
            cells = [Cell.from_str(line) for line in lines]
        return Cells(list(filter(None, cells)))

    def clip_with_bounding_box(self, cells, xmax, ymax, zmax):
        box = (0, xmax, 0, ymax, 0, zmax)
        return Cells([cell for cell in cells if cell.is_inside_box(*box)])

    def clip_enclosed(self, cells, mesh):
        return Cells([cell for cell in cells if cell.is_inside_mesh(mesh)])

    def clip(self, mesh, bounding_box=None):
        """Extract the cells enclosed by either the mesh or bounding box."""
        cells = tqdm(self, desc="Clipping cells")

        if bounding_box:
            return self.clip_with_bounding_box(cells, *bounding_box)
        else:
            return self.clip_enclosed(cells, mesh)

    def hematocrit(self, mesh):
        """Returns an estimated hematocrit value for the given mesh.

        This assumes the approximated red blood cell volume of 90 cubic
        micrometre (µm^3).
        """
        return 90 * len(self.cells) / mesh.volume
