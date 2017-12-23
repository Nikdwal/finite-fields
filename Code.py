import itertools
from util import *
from LinearAlgebra import *
from FiniteField import *

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

class BCHCode:
    def __init__(self, n, field, t, l, beta):
        self.n = n
        self.field = field
        self.t = t
        self.l = l
        self.beta = beta
        self.g = self.make_generator_polynomial()

    def _compute_syndromes(self, v):
        return [None] + [v(self.beta ** (self.l + j - 1)) for j in range(1, 2*self.t + 1)]

    # !! returns Lambda polynomial, not Xk and Yk
    def pgz(self, v, verbose=False):
        beta,l,t = self.beta, self.l, self.t
        S = self._compute_syndromes(v)
        if verbose:
            print("PGZ on v(x) =", v)
        for nu in range(t, 0, -1):
            if verbose:
                print("Trying ν = " + str(nu) + ". ", end="")
            H = Matrix([
                [S[row+col+1] for col in range(nu)]
                for row in range(nu)
            ])
            try:
                Lambda_coef_vector = H.solve(Vector([- S[nu + i + 1] for i in range(nu)]))
                Lambda = Lambda_coef_vector.to_polynomial() + Polynomial.x_to_power(nu, self.field)
                if verbose:
                    print("H is nonsingular.\nResult: Λ(x) =", Lambda, "\n")
                return Lambda
            except ZeroDivisionError:
                # This matrix is singular (the determinant is zero)
                if verbose:
                    print("H is singular.")
                continue
        if verbose:
            print("ν = 0. There were no errors.\n")
        return Polynomial.zero(self.field)

    # verbose output is optional
    def bm(self, v, verbose=False):
        beta = self.beta

        if verbose:
            print("Berlekamp-Massey on v(x) =", str(v), "\n")
            tbl = [["n+d", "Δ", "n", "d", "Λ(x)", "Λ*(x)"],
                   ["==","==", "==", "==", "==", "=="],
                   ["0", "/", "0", "0", "1", "0"]]

        GFq = self.field
        S = [None] + [v(beta ** (self.l + j - 1)) for j in range(1, 2*self.t + 1)]
        Lambda, d, n, Lambda_star = Polynomial.one(GFq), 0, 0, Polynomial.zero(GFq)
        while n + d < 2*self.t:
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
    def forney(self, v, Lambda, verbose=False):
        t,l, beta = self.t, self.l, self.beta
        # determine X_k: Lambda(x) = product(x - X_k)
        X = [root for root in Lambda.find_roots()]

        syndromes = self._compute_syndromes(v)
        S = Polynomial([syndromes[2*t - k] for k in range(0, 2*t)])
        Omega = (S * Lambda) % Polynomial.x_to_power(2*t, Lambda.field)
        Lambdaprime = Lambda.differentiate()

        nu = Lambda.degree()
        Y = [-Omega(X[k]) / (X[k]**(2*t+l) * Lambdaprime(X[k])) for k in range(nu)]
        # if everything went well, Y_k should be in the original field GF(q)

        GFq = self.field
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
        if verbose:
            print("Forney's algorithm on v(x) =", v)
            print("X_k: ", X, "=", ["β^" + str(loc) for loc in error_locations])
            print("S(x) =", S)
            print("Ω(x) =", Omega)
            print("Λ'(x) =", Lambdaprime)
            print("Y_k: = ", Y)
            print("e(x) = ", e, "\n")
        return e

    def make_generator_polynomial(self):
        q = len(self.field)
        beta = self.beta
        cosets = cyclotomic_cosets(q, self.n)
        factors = []
        powers = [i for i in range(self.l, self.l + 2*self.t)]
        x = beta.field.x()
        for coset in cosets:
            # multiply all factors beta**i in this coset if any element i is inbetween l and l+2t-1
            if any([power in coset for power in powers]):
                factors.append(product([x - beta**i for i in coset]))
        g = product(factors)
        return g

    def find_error(self, v, verbose=False):
        Lambda = self.bm(v, verbose)
        e = self.forney(v, Lambda, verbose)
        return e

    def find_error_pgz(self, v, verbose=False):
        Lambda = self.pgz(v, verbose)
        e = self.forney(v, Lambda, verbose)
        return e

    def correct_error(self, v, e, verbose=False):
        c = v - e
        i = c / self.g
        if verbose:
            print("Decoding:")
            print("c(x) = v(x) - e(x) = ", c)
            print("i(x) = c(x) / g(x) = ", i, "\n")
        return i

    def decode(self, v, verbose=False):
        e = self.find_error(v, verbose)
        return self.correct_error(v, e, verbose)

    def decode_pgz(self, v, verbose=False):
        e = self.find_error_pgz(v, verbose)
        return self.correct_error(v, e, verbose)


    def encode(self, i):
        return i * self.g

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
    t,l,v = 3, 2,
    code = BCHCode(n,gf11,t,l,beta)
    Lambda = code.decode(v, verbose=True) # code.berlekamp_massey(v, verbose=True)
    #i = code.decode(v)
    #print(i)

