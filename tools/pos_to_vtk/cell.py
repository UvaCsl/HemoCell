import enum
import pyvista as pv
from scipy.spatial.transform import Rotation

# For the simplified visualisation of both red blood cells and platelets, the
# PyVista ParametricEllipsoids are used, for which the `CELL_RESOLUTION`
# indicates the resolution for the parametric object along each axis. A
# relatively coarse resolution in chosen, to speed up visualisation and
# inside/outside comparisons.
CELL_RESOLUTION = 10


class CellType(enum.Enum):
    """An enumerated type to distinguish cell types."""
    RBC = 1
    PLT = 2

    def geometry(self):
        """Returns the geometry for the given cell type.

        Currently the cell types are all ellipsoidal and are distinguished
        through setting different scaling values in the visualisation. This has
        the additional benefit that the cells can be discriminated easily in
        Paraview based on their given scaling value.
        """
        return pv.ParametricEllipsoid(
            2.2, 4.2, 4.2,      # Human blood
            # 1.7, 2.9, 2.9,    # Mouse blood
            u_res=CELL_RESOLUTION,
            v_res=CELL_RESOLUTION,
            w_res=CELL_RESOLUTION)

    def scaling(self):
        """Returns the scaling value for the cell's geometry.

        This is not an _exact_ scaling, it is an approximate representation
        such that the cells look "OK" on visual inspection.
        """
        if self == CellType.RBC:
            return 1.0

        if self == CellType.PLT:
            return 1./3.

        raise Exception("Scaling is not provided for this cell type.")


class Cell:
    """A simplified representation of a cell.

    The cell only holds the information present in the pos files, i.e. its
    position (center) and rotation angles. The rotation is defined as Euler
    angles with order `x, y, z` as defined by the `packCells` tool.
    """

    def __init__(self, x, y, z, rx, ry, rz):
        self.x = x
        self.y = y
        self.z = z
        self.rx = rx
        self.ry = ry
        self.rz = rz
        self.rotation = Rotation.from_euler('xyz', self.angles(), degrees=True)

    def __repr__(self):
        return f'{self.x} {self.y} {self.z} {self.rx} {self.ry} {self.rz}'

    @staticmethod
    def from_str(string):
        """Parse a `Cell` directly from its string definition."""
        splits = string.split(' ')
        if len(splits) == 6:
            return Cell(*[float(s) for s in splits])

    def center(self):
        """The location of the cell's center."""
        return (self.x, self.y, self.z)

    def angles(self):
        """The Euler rotation angles defined in order `x, y, z`."""
        return (self.rx, self.ry, self.rz)

    def direction(self):
        """Returns a unit vector rotated with the cell's rotation angles.

        The default orientation for the direction vector is assumed aligned
        with the `x` axis, which matches the geometry defined in the `CellType`
        enumerated type.
        """
        return self.rotation.apply([0, 1, 0])

    def bounding_box(self, x=4.2, y=2.2, z=4.2):
        """Returns the bounds of the _rotated_ bounding box of the cell."""
        rotated_box = self.rotation.apply([
            [-x, -y, -z],
            [+x, -y, -z],
            [-x, +y, -z],
            [+x, +y, -z],
            [-x, -y, +z],
            [+x, -y, +z],
            [-x, +y, +z],
            [+x, +y, +z],
        ])

        # Rearrange the coordinates from `[[x0, y0, z0], [x1, y1, z1]], [...]]`
        # to `[(x0, x1, ...), (y0, y1, ...), (z0, z1, ...)]`.
        bounds = list(zip(*rotated_box))

        # Extract the min/max coordinate and shift towards the cell's center.
        return [
            *[c + b for (c, b) in zip(self.center(), map(min, bounds))],
            *[c + b for (c, b) in zip(self.center(), map(max, bounds))],
        ]

    def is_inside_box(self, xmin, xmax, ymin, ymax, zmin, zmax):
        """Returns True when the cell fits fully inside the bounding box."""
        bounding_box = self.bounding_box()
        return (bounding_box[0] > xmin and
                bounding_box[1] > ymin and
                bounding_box[2] > zmin and
                bounding_box[3] < xmax and
                bounding_box[4] < ymax and
                bounding_box[5] < zmax)

    def is_inside_mesh(self, mesh):
        """Returns True when all corners of its bounding box are enclosed."""
        bounding_box = pv.Box(bounds=(self.bounding_box()))
        selected = bounding_box.select_enclosed_points(mesh)
        # The `select_enclosed_points` returns an mask/indicator array with 0/1
        # for points that are outside/inside the mesh respectively. Only when
        # all eight points are inside (i.e. are set to 1) the cell is
        # considered fully inside the mesh.
        return sum(selected['SelectedPoints']) == 8
