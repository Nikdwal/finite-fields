import itertools
from util import max_num_corrected_errors
from LinearAlgebra import Vector
from FiniteField import Polynomial

class LinearBlockCode:
    def __init__(self, G, H=None):
        # G must be a matrix with independent rows
        self.G = G
        self.Gt = G.transpose()
        if H is None:
            H = G.nullspace()
        self.H = H
        self.Ht = H.transpose()
        self.distance = self.Ht.num_independent_rows()
        self._table = self._construct_decode_table()

    def get_field(self):
        return self.generator().get_field()

    def n(self):
        return self.G.width()

    def k(self):
        return self.G.height()

    def _construct_decode_table(self):
        n = self.n()
        F = self.get_field()
        error_positions = range(n)
        table = {}
        nonzeros = [F.generator_to_power(i) for i in range(len(F) - 1)]
        # make every possible error: the parametres that can vary are the number of errors,
        # the positions of those errors, and the bits that make up the errors
        for num_errors in range(max_num_corrected_errors(self.distance) + 1):
            for positions in itertools.combinations(error_positions, num_errors):
                for error_bits in itertools.product(nonzeros, repeat=num_errors):
                    error = [F.zero() for i in range(n)]
                    for i in range(num_errors):
                        error[positions[i]]  = error_bits[i]
                    error = Vector(error)
                    s = self.compute_syndrom(error)
                    table[s] = error
        return table

    def encode(self, word):
        return word * self.G

    def decode_table(self):
        return self._table

    def find_error(self, word):
        return self._table[self.compute_syndrom(word)]

    def generator(self):
        return self.G

    def parity_test_matrix(self):
        return self.H

    def compute_syndrom(self, word):
        return word * self.Ht

# TODO refactor copypasta into a superclass
class CyclicCode:
    def __init__(self, g, n):
        # G must be a matrix with independent rows
        self.g = g
        self.codelength = n
        xnmin1 = Polynomial.x_to_power(n, g.field) - Polynomial.one(g.field)
        # h(x) isn't really used but just store it anyway
        self.h = xnmin1 / g
        self.distance = self.find_distance()
        self._table = self._construct_decode_table()

    def get_field(self):
        return self.generator().field

    def n(self):
        return self.codelength

    def k(self):
        return self.n() - self.generator().degree()

    def find_distance(self):
        d = float("inf")
        for i_coefs in itertools.product(self.get_field(), repeat=self.k()):  # i_coefs is in GF(q) x GF(q) x ...
            i = Polynomial(i_coefs)
            if not i.is_zero():
                c = self.g * i
                c_coefs = Vector(c.coef)
                d = min(d, c_coefs.weight())
        return d

    def _construct_decode_table(self):
        n = self.n()
        F = self.get_field()
        error_positions = range(n)
        table = {}
        nonzeros = [F.generator_to_power(i) for i in range(len(F) - 1)]
        # make every possible error: the parametres that can vary are the number of errors,
        # the positions of those errors, and the bits that make up the errors
        for num_errors in range(max_num_corrected_errors(self.distance) + 1):
            for positions in itertools.combinations(error_positions, num_errors):
                for error_bits in itertools.product(nonzeros, repeat=num_errors):
                    error = [F.zero() for i in range(n)]
                    for i in range(num_errors):
                        error[positions[i]]  = error_bits[i]
                    error = Polynomial(error)
                    s = self.compute_syndrome(error)
                    table[s] = error
        return table

    def encode(self, word):
        return word * self.g

    def decode(self, word):
        v = word
        e = self.find_error(v)
        c = v - e
        i = c / self.g
        return i

    def decode_table(self):
        return self._table

    def find_error(self, word):
        return self._table[self.compute_syndrome(word)]

    def generator(self):
        return self.g

    def parity_test_polynomial(self):
        return self.h

    def compute_syndrome(self, word):
        return word % self.g

if __name__ == "__main__":
    from LinearAlgebra import *
    from FiniteField import *

    GF3 = IntegerField(3)
    x = GF3.x()
    v = 1 + x**2 + 2*x**4
    e = Polynomial(GF3[0,1,0,0,2,0,0,0])
    g = (2 + x + x ** 2) * (1 + x ** 2) * (1 + x)
    n = 8
    code = CyclicCode(g, 8)
