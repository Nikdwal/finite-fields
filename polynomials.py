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
            # p(x) * (a + bx + cx^2) = a * p(x) + bx * p(x) + cx^2 *p(x)
            prod = Polynomial.zero(self.field)
            dq = q.degree()
            for i in range(len(other)):
                prod += q[i] * p.multiply_monic_mononomial(i)
            return prod

    # multiply by a scalar
    def scale(self, scalar):
        assert type(scalar) is type(self.coef[0])
        coefs = [scalar * c for c in self.coef]
        return Polynomial(coefs, self.field)

    def __divmod__(self, other):
        assert self.field == other.field
        return divide_polynomials(self, other)

    def __floordiv__(self, other):
        assert self.field == other.field
        return divmod(self, other)[0]

    def __mod__(self, other):
        assert self.field == other.field
        return divmod(other, self)[1]

    def __truediv__(self, other):
        assert self.field == other.field
        quotient, remainder = divmod(self, other)
        if remainder.is_zero():
            return quotient
        else:
            raise ValueError("Cannot divide these two polynomials. The remainder is non-zero.")

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
print(p)
print(q)
print(two * p)
print(p * q)
print((p * q) / p)
print(p % ((p * q) + Polynomial.monic_mononomial(1, Z5)))
print(p - q)