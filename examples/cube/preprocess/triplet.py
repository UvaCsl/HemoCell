from utilities import even_repartition


class Triplet:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __iter__(self):
        for v in (self.x, self.y, self.z):
            yield v

    def __add__(self, other):
        return Triplet(self.x + other.x, self.y + other.y, self.z + other.z)

    def __mul__(self, other):
        if isinstance(other, float) or isinstance(other, int):
            return Triplet(self.x * other, self.y * other, self.z * other)
        return Triplet(self.x * other.x, self.y * other.y, self.z * other.z)

    def __repr__(self):
        return f'{self.x}x{self.y}x{self.z}'


class AtomicBlock(Triplet):
    ...


class Pattern(Triplet):
    ...


class Decomposition(Triplet):
    def __init__(self, n_cpu):
        super().__init__(*even_repartition(n_cpu, 3))
