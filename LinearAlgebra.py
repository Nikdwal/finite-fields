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
        if length <= l:
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
        return self.get_elems().__iter__()

    def __getitem__(self, m):
        # use the fact that if a polynomial of length l < L, p[L] returns zero
        return self.polynomial[m]

    def __repr__(self):
        return repr(self.get_elems())

    def __eq__(self, other):
        return isinstance(other, Vector) and self.polynomial == other.polynomial and len(self) == len(other)

    def __hash__(self):
        # TODO: is this correct?
        return hash(self._attrs)

    def __neg__(self):
        return Vector.from_polynomial(- self.polynomial, len(self))

    def scale(self, scalar):
        # you can only multiply a vector by a scalar
        assert isinstance(scalar, FieldElement) and scalar.field == self.get_field()

        return Vector.from_polynomial(scalar * self.polynomial, len(self))

    def div(self, scalar):
        assert isinstance(scalar, FieldElement) and scalar.field == self.get_field()
        return self.scale(scalar ** -1)

    def dot_product(self, other):
        assert isinstance(other, Vector)
        return Vector(dot_product(self.get_elems(), other.get_elems()))

    def __add__(self, other):
        if not isinstance(other, Vector):
            raise TypeError("Vectors can only be added to vectors.")
        if len(other) != len(self):
            raise ValueError("Cannot add vectors of different lengths")
        return Vector.from_polynomial(self.polynomial + other.polynomial, len(self))

    def __sub__(self, other):
        return self + (-other)

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

    # "rows" must be a list of lists of FieldElems
    def __init__(self, rows):
        width = len(rows[0])
        field = rows[0][0].field
        assert all([len(row) == width for row in rows])
        assert all([isinstance(e, FieldElement) and e.field == field for e in itertools.chain.from_iterable(rows)])
        self.rows = rows

    def cast_to_field(self, field):
        return Matrix([[field.cast(elem) for elem in row] for row in self])

    def is_zero(self):
        return all([all([elem.is_zero() for elem in row]) for row in self])

    def __eq__(self, other):
        return (self - other).is_zero()

    def __hash__(self):
        return hash(self.rows())

    def __eq__(self, other):
        return

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
        assert isinstance(other, Matrix)
        assert(self.width() == other.height())
        F = self.get_field().largest(other.get_field())
        m1 = self.cast_to_field(F)
        m2 = other.cast_to_field(F)

        height = m1.height()
        width = m2.width()
        m2_transpose = m2.transpose()
        # The ij'th element is the dot product of the i'th row in m1 and the j'th row in m2^T
        prod = [[None for j in range(width)] for i in range(height)]
        for i in range(height):
            for j in range(width):
                prod[i][j] = dot_product(m1[i], m2_transpose[j])
        return Matrix(prod)

    def __rmul__(self, other):
        if isinstance(other, FieldElement) or isinstance(other, int):
            return Matrix([[other * elem for elem in row] for row in self])
        elif isinstance(other, Vector):
            row_matrix = Matrix.from_vector(other)
            product = row_matrix * self
            return Vector(product[0])
        else:
            raise ValueError("Cannot multiply these items.")

if __name__ == "__main__":
    from FiniteField import IntegerField
    Z2 = IntegerField(2)
    z, o = Z2.zero(), Z2.one()
    m = Matrix([[o,z], [z,o]])
    v = Vector([z,o])
    print(v, "\n")
    print(m, "\n")
    print(v*m)
    pass