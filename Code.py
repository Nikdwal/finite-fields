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

# === BCH codes

# verbose output is optional
def berlekamp_massey(t, l, v, beta, verbose=False):
    if verbose:
        print("Berlekamp-Massey on v(x) =", str(v), ":\n")
        tbl = [["n+d", "Δ", "n", "d", "Λ(x)", "Λ*(x)"],
               ["==","==", "==", "==", "==", "=="],
               ["0", "/", "0", "0", "1", "0"]]

    GFq = v.field
    S = [None] + [v(beta ** (l + j - 1)) for j in range(1, 2*t + 1)]
    Lambda, d, n, Lambda_star = Polynomial.one(GFq), 0, 0, Polynomial.zero(GFq)
    while n + d < 2*t:
        Delta = sum([Lambda[m] * S[m+n+1] for m in range(0, d+1)])
        if Delta.is_zero():
            n += 1
        elif n < d:
            Lambda -= Delta * Lambda_star.multiply_x_to_power(d - n - 1)
            n += 1
        else:
            (Lambda, Lambda_star) = (Lambda.multiply_x_to_power(n-d+1) - Delta * Lambda_star, Lambda / Delta)
            (n, d) = (d, n + 1)
        if verbose:
            tbl.append([str(x) for x in [n+d, Delta, n, d, Lambda, Lambda_star]])

    if verbose:
        tbl = [[entry + " | " for entry in row] for row in tbl]
        column_widths = [max(len(item) for item in [row[i] for row in tbl]) for i in range(len(tbl[0]))]
        format = ''.join(['%%%ds' % width for width in column_widths])
        for row in tbl:
            print(format % tuple(row))
        print("\nResult: Λ(x) = ", Lambda, "\n")

    return Lambda

# return the error
def forney(t, l, v, beta, Lambda):
    # determine X_k: Lambda(x) = product(x - X_k)
    X = [root for root in Lambda.find_roots()]

    S_coefs = [v(beta ** (l + j - 1)) for j in range(2 * t, 0, -1)]
    S = Polynomial(S_coefs)
    Omega = (S * Lambda) % Polynomial.x_to_power(2*t, Lambda.field)
    Lambdaprime = Lambda.differentiate()
    nu = Lambda.degree()
    Y = [-Omega(X[k]) / (X[k]**(2*t+l) * Lambdaprime(X[k])) for k in range(nu)]
    # if everything went well, Y_k should be in the original field GF(q)

    GFq = v.field
    F = beta.field

    # WARNING: the following code block is very complicated and distracting
    # Now we want to determine the error locations i_k as X[k] = beta^i_k
    # The problem is that we do not have a lookup table for powers of beta, but we do have a table
    # for powers of alpha (in the extended field F)
    # say beta = alpha^r. And X[k] = alpha^s = beta^i_k. This means alpha^s = alpha^(r*i_k)
    N = len(F) - 1
    # Therefore s = r*i_k (mod N) because <F\{0},*> is a cyclic group with generator alpha
    # In other words, s = r*i_k + j*N for some (positive or negative) integer j
    # this means i_k = (s - j*N)/r
    # this reduces the search for i_k to just a few values for j (which is expected to be quite small: typically j=0)
    # Depending on whether you already have a lookup table for powers of alpha but no powers of beta, it might
    # not be a good idea to use this method when computing this by hand
    error_locations = [None for x in X]
    for k in range(len(X)):
        s = F.generator_exponent(X[k])
        r = F.generator_exponent(beta)
        j = 0
        ik = 0
        while beta ** ik != X[k]:
            ik = (s - j*N) // r
            # sequence: 0,1,-1,2,-2,3,-3,...
            if j > 0:
                j = -j
            else:
                j = -j + 1

        error_locations[k] = ik

    e = sum(Y[k] * Polynomial.x_to_power(error_locations[k], GFq) for k in range(nu))
    return e

# === get_generator_polynomial
def generator_polynomial(n,q,t,l,beta):
    cosets = cyclotomic_cosets(q, n)
    factors = []
    powers = [i for i in range(l, l + 2*t)]
    x = beta.field.x()
    for coset in cosets:
        # multiply all factors beta**i in this coset if any element i is inbetween l and l+2t-1
        if any([power in coset for power in powers]):
            factors.append(product([x - beta**i for i in coset]))
    g = product(factors)
    return g

if __name__ == "__main__":
    from LinearAlgebra import *
    from FiniteField import *

    # GF3 = IntegerField(3)
    # x = GF3.x()
    # v = 1 + x**2 + 2*x**4
    # e = Polynomial(GF3[0,1,0,0,2,0,0,0])
    # g = (2 + x + x ** 2) * (1 + x ** 2) * (1 + x)
    # n = 8
    # code = CyclicCode(g, 8)


    q = 11
    gf11 = IntegerField(q)
    w = Polynomial(gf11[7,1,1])
    n = 15
    gf121 = ExtendedField(gf11,2,"α", w)
    beta = gf121.generator_to_power(8)
    t,l,v = 3, 2, Polynomial(gf11[3, 5, 4, 3, 6, 9, 10, 9, 8, 10, 8, 4, 6, 4, 7])
    Lambda = berlekamp_massey(t,l,v, beta, verbose=True)
    e = forney(t,l,v,beta,Lambda)
    x = gf121.x()
    g = generator_polynomial(n,q,t,l,beta)
    c = v - e
    i = c / g
    print(i)

