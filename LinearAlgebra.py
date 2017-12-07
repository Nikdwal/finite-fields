from FiniteField import Polynomial, FieldElement
from util import dot_product
import itertools

# A **row** vector
class Vector():
    def __init__(self, components):
        # Store each vector as a polynomial.
        # Complications arise when higher degree terms are zero. The constructor
        # for polynomials would throw these away because these are inconvenient to work with.
        # Since vectors must have a fixed length, we must store this as a variable in addition to the polynomial.
        self.polynomial = Polynomial(components)
        self._length = len(components)
        self._attrs = (self.polynomial, self._length)

    @staticmethod
    def from_polynomial(polynomial, length):
        l = len(polynomial)
        if length <= polynomial:
            return Vector(polynomial.coef[:length])
        else:
            zero = polynomial.field.zero()
            return Vector(polynomial.coef + [zero for i in range(length - l)])

    def __len__(self):
        return self._length

    def get_field(self):
        return self.polynomial.field

    def trailingzeros(self):
        # zeros at the end of the vector that are omitted by the polynomial representation
        return [self.get_field().zero() for i in range(len(self) - len(self.polynomial))]

    def get_elems(self):
        return self.polynomial.coef + self.trailingzeros()

    def __iter__(self):
        self.get_elems().__iter__()

    def __getitem__(self, m):
        # use the fact that if a polynomial of length l < L, p[L] returns zero
        return self.polynomial[m]

    def __repr__(self):
        return repr(self.get_elems())

    def __eq__(self, other):
        return self.polynomial == other.polynomial and len(self) == len(other)

    def __hash__(self):
        # TODO: is this correct?
        return hash(self._attrs)

    def __neg__(self):
        return Vector.polynomial(- self.polynomial, len(self))

    def __rmul__(self, scalar):
        # you can only multiply a vector by a scalar
        assert type(scalar) is FieldElement and scalar.field == self.get_field()

        return Vector.from_polynomial(scalar * self.polynomial, len(self))

    def dot_product(self, other):
        assert type(other) is Vector
        return Vector(dot_product(self.get_elems(), other.get_elems()))

    # TODO: make it possible to divide by a scalar by first implementing this functionality for polynomials

    def __add__(self, other):
        assert type(other) is Vector and other.get_field() == self.get_field() and len(other) == len(self)
        return Vector.from_polynomial(self.polynomial + other.polynomial, len(self))

    def __sub__(self, other):
        assert type(other) is Vector and other.get_field() == self.get_field() and len(other) == len(self)
        return Vector.from_polynomial(self.polynomial - other.polynomial, len(self))

    def weight(self):
        return len([c for c in self.polynomial if not c.is_zero()])

    def distance(self, other):
        return (self - other).weight()

    def __gt__(self, other):
        return self.weight() > other.weight

    def __lt__(self, other):
        return self.weight() < other.weight

    def __le__(self, other):
        return self.weight() <= other.weight

    def __ge__(self, other):
        return self.weight() >= other.weight

class Matrix:
    # rows must be a list of lists of FieldElems
    def __init__(self, rows):
        width = len(rows[0])
        field = rows[0][0].field
        assert all([len(row) == width for row in rows])
        assert all([type(e) is FieldElement and e.field == field for e in itertools.chain.from_iterable(rows)])
        self.rows = rows

    # make a **row** vector
    @staticmethod
    def from_vector(vector):
        return Matrix([vector.get_elems()])

    def __repr__(self):
        s = ""
        for row in self:
            s += repr(row) + "\n"
        return s

    def get_field(self):
        return self.rows[0][0].field

    def height(self):
        return len(self.rows)

    def width(self):
        return len(self.rows[0])

    # get the row_number'th row
    def __getitem__(self, row_number):
        return self.rows[row_number]

    # get an iterator of the rows
    def __iter__(self):
        return self.rows.__iter__()

    def transpose(self):
        trans = [[None for j in range(self.height())] for i in range(self.width())]
        for c in range(self.width()):
            for r in range(self.height()):
                trans[c][r] = self[r][c]
        return Matrix(trans)

    def __add__(self, other):
        assert self.height() == other.height() and self.width() == other.width()
        return Matrix([(Vector(self[r]) + Vector(other[r])).get_elems() for r in range(self.height())])

    def __neg__(self):
        return Matrix([(- Vector(row)).get_elems() for row in self])

    def __sub__(self, other):
        return self + (- other)

    def __mul__(self, other):
        if type(other) is Matrix:
            assert(self.width() == other.height() and self.get_field() == other.get_field())
            height = self.height()
            width = other.width()
            other_transpose = other.transpose()
            # ij'th element is the dot product of the i'th row in self and the j'th row in other^T
            prod = [[None for j in range(width)] for i in range(height)]
            for i in range(height):
                for j in range(width):
                    prod[i][j] = dot_product(self[i], other_transpose[j])
            return Matrix(prod)

    def __rmul__(self, other):
        if type(other) is FieldElement:
            assert other.field == self.get_field()
            return Matrix([(other * Vector(row)).get_elems() for row in self])

        elif type(other) is Vector:
            row_matrix = Matrix.from_vector(other)
            product = row_matrix * self
            return Vector(product[0])

if __name__ == "__main__":
    # from FiniteField import IntegerField
    # Z2 = IntegerField(2)
    # z, o = Z2.zero(), Z2.one()
    # m = Matrix([[o,z], [z,o]])
    # v = Vector([z,o])
    # print(v, "\n")
    # print(m, "\n")
    # print(v*m)
    pass