class Cell:
    def __init__(self, x, y, z, rx, ry, rz):
        """Basic representation of a cell using its position and rotation."""
        self.x = x
        self.y = y
        self.z = z
        self.rx = rx
        self.ry = ry
        self.rz = rz

    def __repr__(self):
        return f'{self.x} {self.y} {self.z} {self.rx} {self.ry} {self.rz}'

    def translate(self, x, y, z):
        return Cell(self.x + x, self.y + y, self.z + z,
                    self.rx, self.ry, self.rz)

    @staticmethod
    def from_str(string):
        splits = string.split(' ')
        assert len(splits) == 6, "Needs 6 values per Cell."
        return Cell(*(float(s) for s in splits))
