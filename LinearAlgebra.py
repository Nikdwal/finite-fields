from FiniteField import Polynomial, FieldElement
from util import *
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

    def to_polynomial(self):
        return self.polynomial

    def __len__(self):
        return self._length

    def get_field(self):
        return self.polynomial.field

    def is_zero(self):
        return self.polynomial.is_zero()

    def cast_to_field(self, field):
        return Vector([field.cast(c) for c in self])

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
        return Vector.from_polynomial(scalar * self.polynomial, len(self))

    def __rmul__(self, other):
        if isinstance(other, int) or isinstance(other, FieldElement):
            return self.scale(other)

    def div(self, scalar):
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
    # value_semantics controls whether the array should be copied (recommended)
    def __init__(self, rows, value_semantics=True):
        width = len(rows[0])
        field = rows[0][0].field
        assert all([len(row) == width for row in rows])
        assert all([isinstance(e, FieldElement) and e.field == field for e in itertools.chain.from_iterable(rows)])
        if value_semantics:
            self.rows = [list(row) for row in rows]
        else:
            # don't bother copying the rows
            self.rows = rows

    def cast_to_field(self, field):
        return Matrix([[field.cast(elem) for elem in row] for row in self])

    def is_zero(self):
        return all([all([elem.is_zero() for elem in row]) for row in self])

    def __eq__(self, other):
        return (self - other).is_zero()

    def __hash__(self):
        return hash(self.rows())

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

    # !! be careful with this
    # It might break integrity
    def __setitem__(self, row_index, value):
        self.rows[row_index] = value

    # get an iterator of the rows
    def __iter__(self):
        return self.rows.__iter__()

    def transpose(self):
        trans = [[None for j in range(self.height())] for i in range(self.width())]
        for c in range(self.width()):
            for r in range(self.height()):
                trans[c][r] = self[r][c]
        return Matrix(trans, value_semantics=False)

    def __add__(self, other):
        assert self.height() == other.height() and self.width() == other.width()
        return Matrix([add_arrays(self[r], other[r]) for r in range(self.height())])

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
        prod = [[None for j in range(width)] for i in range(height)]
        for i in range(height):
            for j in range(width):
                prod[i][j] = dot_product(m1[i], [m2[x][j] for x in range(m2.height())])
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

    def is_square(self):
        return self.width() == self.height()

    def rowswap(self, row1, row2):
        tmp = self[row1]
        self[row1] = self[row2]
        self[row2] = tmp

    def rowmul(self, rownumber, scalar):
        for i in range(self.width()):
            self[rownumber][i] = scalar * self[rownumber][i]

    def add_multiple_of_row(self, row1, row2, scalar):
        for i in range(self.width()):
            self[row1][i] = self[row1][i] + scalar * self[row2][i]

    def to_echelon_form(self):
        U = self
        height = U.height()
        width = U.width()

        row = 0
        for col in range(width):
            if row >= height:
                break
            pivot = U[row][col]
            if pivot.is_zero():
                found_nonzero = False
                for j in range(row+1, height):
                    if not U[j][col].is_zero():
                        # swap rows
                        U.rowswap(row, j)
                        found_nonzero = True
                        pivot = U[row][col]
                        break
                if not found_nonzero:
                    # don't update the row index, only update the column index
                    continue
            for j in range(row+1, height):
                U.add_multiple_of_row(j, row, - U[j][col] / pivot)
            row += 1

    def to_reduced_echelon_form(self):
        U = self
        height = U.height()
        width = U.width()

        row = 0
        for col in range(width):
            if row >= height:
                break
            pivot = U[row][col]
            if pivot.is_zero():
                found_nonzero = False
                for j in range(row + 1, height):
                    if not U[j][col].is_zero():
                        # swap rows
                        U.rowswap(row, j)
                        found_nonzero = True
                        pivot = U[row][col]
                        break
                if not found_nonzero:
                    # don't update the row index, only update the column index
                    continue
            # eliminate
            U.rowmul(row, pivot ** -1)
            for j in range(height):
                if j != row:
                    U.add_multiple_of_row(j, row, - U[j][col])
            row += 1

    def echelon_form(self):
        U = Matrix(self.rows)
        U.to_echelon_form()
        return U

    def reduced_echelon_form(self):
        U = Matrix(self.rows)
        U.to_reduced_echelon_form()
        return U

    def product_diagonal_elements(self):
        product = self.get_field().zero()
        for i in range(min(self.height(), self.width())):
            product += self[i][i]
        return product

    def rank(self):
        U = self.echelon_form()
        rank = 0
        for row in U:
            if all([u.is_zero() for u in row]):
                break
            rank += 1
        return rank

    def nullity(self):
        return self.width() - self.rank()

    # solve Ax = b with b a (column) vector
    # or equivalently: xA^T = b with b a row vector
    def solve(self, righthandside):
        assert isinstance(righthandside, Vector)
        if not self.is_square():
            raise NotImplementedError("Only square systems have been implemented so far")

        if len(righthandside) != self.height():
            raise ValueError("Cannot solve this system. The dimensions don't match.")
        A = Matrix(self.rows)
        for row in range(A.height()):
            A[row].append(righthandside[row])
        A.to_echelon_form()

        if A.product_diagonal_elements().is_zero():
            raise ZeroDivisionError("Error. This matrix is singular.")

        # back substitution
        # remember that b is the last column of A
        n = A.height()
        F = A.get_field()
        x = [F.zero() for i in range(n)]
        x[-1] = 1 / A[-1][-2] * A[-1][-1]
        for i in range(n - 2, -1, -1):
            x[i] = (A[i][-1] - dot_product(A[i][i+1 : n], x[i+1:])) / A[i][i]
        return Vector(x)

    # Return a matrix whose row space is the null space of this matrix.
    # In other words, if self is G and nullspace is H, then G * H^T = 0
    def nullspace(self):
        # This boils down to solving a homogeneous system
        U = self.echelon_form()
        if U.width() < U.height():
            # This probably won't be too difficult, but you have to think about all the edge cases.
            raise NotImplementedError("This matrix shape is not supported yet.")

        F = U.get_field()
        if U.is_square() and not U.product_diagonal_elements().is_zero():
            # If the matrix is nonsingular, the null space is just zero.
            return Matrix([[F.zero() for i in range(U.width())]], value_semantics=False)

        # Backward substitution
        row = U.height() - 1
        H = [[F.zero() for i in range(U.width())]]
        H[0][-1] = F.one()
        for col in range(U.width()-2, -1, -1):
            if row < 0 or any(not U[row][j].is_zero() for j in range(col)):
                # add one dimension to the null space by adding a copy of the last vector
                # where the col'th component is an arbitrary nonzero number
                # all other vectors will have a zero here
                H.append(H[-1].copy())
                H[-1][col] = F.one()
            else:
                for vector in H:
                    vector[col] = - dot_product(U[row][col+1:], vector[col+1:]) / U[row][col]
                row -= 1

        return Matrix(H, value_semantics=False)

    # THE MATRIX MUST BE IN THE FORM [ I | P] BEFORE YOU USE THIS
    def nullspace_from_standard_form(self):
        # assuming the matrix is in a form [ I (n x n) | P (n x (m-n)) ]
        # the null space is given by [- P (n x (m-n)) ^T | I (m-n) x (m-n))]
        n = self.height()
        m = self.width()
        F = self.get_field()
        H = [
            [- self[j][n + i] for j in range(n)] # - P^T
            + [F.zero() for k in range(i)] + [F.one()] + [F.zero() for k in range(m-n-i-1)] # identity
            for i in range(m - n) # i is the row index of the result
        ]
        return Matrix(H, value_semantics=False)

    def is_singular(self):
        return not self.is_square() or self.echelon_form().product_diagonal_elements().is_zero()

    # Return the greatest number d such any combinations of d or fewer rows is linearly independent
    def num_independent_rows(self):
        # TODO: is there a cleverer way than this brute force approach?
        # e.g. doing row-ops first will generally not work (it would change the answer)

        rank_upperbound = min(self.height(), self.width())
        for n in range(2, rank_upperbound):
            for row_collection in itertools.combinations(self.rows, n):
                S = Matrix(row_collection, value_semantics=False)
                # rank computation will not modify the matrix so we don't have to copy the rows in the constructor
                if S.rank() < n:
                    return n - 1
        return rank_upperbound

if __name__ == "__main__":
    from FiniteField import IntegerField
    Z2 = IntegerField(2)
    H = Matrix(Z2[[1,0,1,0,1,0,1],[0,1,1,0,0,1,1],[0,0,0,1,1,1,1]])
    Ht = H.transpose()
    print(max_num_corrected_errors(Ht.num_independent_rows()))
    pass