import math
import abc

from math import prod
from triplet import Pattern
from utilities import cubic_root, doubling_range
from utilities import is_exact_cube, even_repartition


class Layout(abc.ABC):
    def __init__(self, min_node=1, max_node=math.inf):
        self.min_node = min_node
        self.max_node = max_node

    @abc.abstractmethod
    def __call__(self, x):
        """Return the current layout in for size `x` in Triplet form."""
        ...

    def __iter__(self):
        for x in doubling_range(start=self.min_node):
            pattern = self(x)
            if prod(pattern) > self.max_node:
                break
            yield pattern


class Linear(Layout):
    def __call__(self, x):
        """Return a linear pattern along the x-axis."""
        return Pattern(x, 1, 1)


class Cubic(Layout):
    def __init__(self, start=1, *args, **kwargs):
        assert is_exact_cube(start) > 0, "Needs to start with a cubic root"
        super().__init__(cubic_root(start), *args, **kwargs)

    def __call__(self, x):
        """Return a cubic pattern in x, y, z."""
        return Pattern(x, x, x)


class Box(Layout):
    def __init__(self, dim, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.previous = None
        self.count = 0
        self.dim = dim

    def __call__(self, x):
        """Return increasing domain along alternating dimensions."""
        if self.previous is None:
            self.previous = Pattern(*even_repartition(x, 3))
            return self.previous

        if self.count == 0:
            self.previous.y *= 2
        elif self.count == 1:
            self.previous.z *= 2
        else:
            self.previous.x *= 2

        self.count = (self.count + 1) % self.dim
        return self.previous


class Planar(Box):
    def __init__(self, *args, **kwargs):
        """A two-dimensional alternating pattern."""
        super().__init__(2, *args, **kwargs)


class Hexahedral(Box):
    def __init__(self, *args, **kwargs):
        """A three-dimensional alternating pattern."""
        super().__init__(3, *args, **kwargs)
