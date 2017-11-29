class Polynomial:
    # coef is a list of coefficients from lowest to highest degree terms
    def __init__(self, coef, field):
        self.field = field

        # remove higher degree terms if their coefficients are zero (those would be redundant)
        first_zero = len(coef)
        for i in range(len(coef) -1, -1, -1):
            if coef[i].is_zero():
                first_zero -= 1
            else:
                break
        self.coef = coef[:first_zero]

    @staticmethod
    def zero(field):
        return Polynomial([], field)

    @staticmethod
    def one(field):
        return Polynomial.monic_mononomial(0, field)

    @staticmethod
    def monic_mononomial(degree, field):
        return Polynomial([field.zero() for i in range(degree)] + [field.one()], field)

    def __str__(self):
        if len(self) == 0:
            return str(self.field.zero())

        superscript = str.maketrans("0123456789", "⁰ ²³⁴⁵⁶⁷⁸⁹")
        s = "" if self[0].is_zero() else str(self[0])
        for i in range(1, len(self)):
            if not self.coef[i].is_zero():
                s += ("" if s == "" else " + ") +\
                     ( "" if self.coef[i].is_one() else str(self.coef[i])) +\
                     "x" + str(i).translate(superscript)
        return s

    # compute the value of the polynomial when substituting x=args[0]
    def __call__(self, *args, **kwargs):
        x = args[0]
        res = self[-1]
        for i in range(self.degree() - 1, -1, -1):
            res = self[i] + res * x
        return res

    # get the n'th degree coefficient
    # usage: p[n] where p is a polynomial
    def __getitem__(self, m):
        return self.coef[m]

    def __eq__(self, other):
        return self.field == other.field and (self - other).is_zero()

    def __hash__(self):
        # TODO: is this correct? does a == b => hash(a) == hash(b) in finite fields?
        return hash(tuple(self.coef))

    def is_zero(self):
        return self.coef == []

    def is_one(self):
        return len(self) == 1 and self.coef[0].is_one()

    def __len__(self):
        return len(self.coef)

    def degree(self):
        return len(self) - 1

    def __add__(self, other):
        assert self.field == other.field
        d1 = self.degree()
        d2 = other.degree()
        if d1 > d2:
            c = other.coef + [self.field.zero() for i in range(d1 - d2)]
            return Polynomial(add_arrays(self.coef, c), self.field)
        else: 
            c = self.coef + [self.field.zero() for i in range(d2 - d1)]
            return Polynomial(add_arrays(other.coef, c), self.field)

    def __neg__(self):
        return Polynomial([-c for c in self.coef], self.field)

    def __sub__(self, other):
        assert self.field == other.field
        return self + (-other)

    # multiply this polynomial by (x ^ power)
    def multiply_monic_mononomial(self, degree):
        return Polynomial([self.field.zero() for i in range(degree)] + self.coef, self.field)

    def __mul__(self, other):
        if self.is_zero():
            return self
        elif type(other) is type(self.coef[0]):
            # in other words: multiply this polynomial by a scalar that is in the same field as the coefficients
            return self.scale(other)
        else:
            assert type(other) is Polynomial and other.field == self.field
            p, q = self, other
            # p(x) * (q_0 + q_1 x + q_2 x^2) = q_0 * p(x) + q_1 x * p(x) + q_2 x^2 * p(x)
            prod = Polynomial.zero(self.field)
            for i in range(len(other)):
                prod += q[i] * p.multiply_monic_mononomial(i)
            return prod

    # multiply by a scalar
    def scale(self, scalar):
        assert type(scalar) is type(self.coef[0])
        coefs = [scalar * c for c in self.coef]
        return Polynomial(coefs, self.field)

    # return the quotient and the remainder
    def __divmod__(self, other):
        assert self.field == other.field
        return divide_polynomials(self, other)

    # return the quotient, but not the remainder
    def __floordiv__(self, other):
        assert self.field == other.field
        return divmod(self, other)[0]

    # return the remainder for other / self (== self % other)
    def __mod__(self, other):
        assert self.field == other.field
        return divmod(other, self)[1]

    # return the quotient iff the remainder is zero
    def __truediv__(self, other):
        # TODO: make it possible to divide by a scalar
        assert self.field == other.field
        quotient, remainder = divmod(self, other)
        if remainder.is_zero():
            return quotient
        else:
            raise ValueError("Cannot divide these two polynomials. The remainder is non-zero.")

    # find an x such that p(x) is zero
    def find_root(self):
        # there is no clever way to do this in finite fields
        for x in self.field:
            if self(x).is_zero():
                return x
        return None

    # @staticmethod
    # def all_monic_polynomials(degree, field):
    #     if degree == 0:
    #         yield Polynomial([field.one()], field)
    #     else:
    #         for c in field:
    #             for p in Polynomial.all_monic_polynomials(degree - 1, field):
    #                 yield p.multiply_monic_mononomial(1) + Polynomial([c], field)

    @staticmethod
    def _make_larger_irreducible_polyns(field, irreducibles):
        # there are a few relatively clever tricks to find all of these given a certain field
        # however, none of them are quite remarkable; nor are they easy to implement algorithmically
        for p in irreducibles:
            for c in field:
                # this will make all monic polynomials of degree one higher
                w = p.multiply_monic_mononomial(1) + Polynomial([c], field)

                if not w.is_reducible():
                    yield w

    def factor(self):
        # TODO: this part is very messy and inefficient
        F = self.field

        a_n = self[-1]
        if a_n.is_one():
            factors = []
            w = self
        else:
            # split off the leading coefficient
            factors = [Polynomial([a_n], F)]
            w = (F.one() / a_n) * self

        irreducibles = [Polynomial([F.one()], F)]
        for deg in range(1, 1 + self.degree() // 2):
            for p in Polynomial._make_larger_irreducible_polyns(F, irreducibles):

                if w.is_one():
                    # completely factorised
                    return factors

                quot, rem = divmod(w, p)
                if rem.is_zero():
                    # found a factor
                    factors.append(p)
                    w = quot

        # TODO: will this even occur?
        return factors + [w]


    def is_reducible(self):
        # TODO: remove copypasta

        F = self.field

        a_n = self[-1]
        if a_n.is_one():
            factors = []
            w = self
        else:
            # split off the leading coefficient
            factors = [a_n]
            w = (F.one() / a_n) * self

        irreducibles = [Polynomial([F.one()], F)]
        for deg in range(1, 1 + self.degree() // 2):
            for p in Polynomial._make_larger_irreducible_polyns(F, irreducibles):
                quot, rem = divmod(w, p)
                if rem.is_zero():
                    # found a factor
                    return True
        return False



def add_arrays(lst1, lst2):
    assert len(lst1) == len(lst2)
    s = [x for x in lst1]
    for i in range(len(lst1)):
        s[i] += lst2[i]
    return s

# returns the quotient and the remainder
def divide_polynomials(numerator, denominator):
    assert numerator.field == denominator.field
    F = numerator.field
    d1 = numerator.degree()
    d2 = denominator.degree()
    if d1 < d2:
        return Polynomial.zero(F), numerator
    deg_quotient = d1 - d2

    # long division ("staartdeling")
    leading_coef = numerator[d1] / denominator[d2]
    quotient_leading_terms = leading_coef * Polynomial.monic_mononomial(deg_quotient, F)
    denom_times_quot = (leading_coef * denominator).multiply_monic_mononomial(deg_quotient)
    smaller_polynomial = numerator - denom_times_quot # The coefficient on the leading term will disappear
    if smaller_polynomial.is_zero():
        return quotient_leading_terms, Polynomial.zero(F)
    quotient_smaller_polynomial, remainder = divide_polynomials(smaller_polynomial, denominator)
    return quotient_leading_terms + quotient_smaller_polynomial, remainder


#===== testing
from FiniteField import PrimeField
Z5 = PrimeField(5)
zero = Z5.zero()
one = Z5.one()
two = one + one
three = one + two
four = one + three
p = Polynomial([one, three, two], Z5)
q = Polynomial([two, four], Z5)
# print(p)
# print(q)
# print(two * p)
# print(p * q)
# print((p * q) / p)
# print(p % ((p * q) + Polynomial.monic_mononomial(1, Z5)))
# print(p - q)
# print(p(four))
#
# factors = p.factor()
# for factor in factors:
#     print(factor)
#
# pol = Polynomial.monic_mononomial(0, Z5)
# for factor in factors:
#     pol *= factor
#
# print(pol)
