import itertools
from util import max_num_corrected_errors
from LinearAlgebra import Vector

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
        zerovector = Vector([F.zero() for i in range(n)])
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
                    s = error * self.Ht
                    table[s] = error
        return table

    def encode(self, word):
        return word * self.G

    def decode_table(self):
        return self._table

    def find_error(self, word):
        return self._table[word * self.Ht]

    def generator(self):
        return self.G

    def parity_test_matrix(self):
        return self.H

if __name__ == "__main__":
    from LinearAlgebra import *
    from FiniteField import *

    GF3 = IntegerField(3)
    G = Matrix(GF3[[2, 0, 1, 1, 2, 1, 0, 0], [0, 2, 0, 1, 1, 2, 1, 0], [0, 0, 2, 0, 1, 1, 2, 1]])
    code = LinearBlockCode(G)
    word = Vector(GF3[1,0,0])
    print(code.decode_table())
