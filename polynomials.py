class Polynomial:
    # coef is a list of coefficients from highest to lowest degree terms
    def __init__(self, coef, field):
        self.field = field
        first_nonzero = -1
        for i in range(len(coef)):
            if not coef[i].is_zero():
                first_nonzero = i
                break
        self.coef = coef[first_nonzero:] if first_nonzero >= 0 else []

    @staticmethod
    def zero(field):
        return Polynomial([], field)

    @staticmethod
    def monic_mononomial(degree, field):
        return Polynomial([field.one()] + [field.zero() for i in range(degree)], field)

    def __str__(self):
        if len(self) == 0:
            return "0"
        s = ""
        d = self.degree()
        sup = str.maketrans("0123456789", "⁰ ²³⁴⁵⁶⁷⁸⁹")
        for i in range(d):
            if not self.coef[i].is_zero():
                s += (str(self.coef[i]) if self.coef[i] != self.field.one() else "") + "x" + str(d - i).translate(sup) + " + "
        s += str(self.coef[d])
        return s

    # get the n'th degree coefficient
    # usage: p[n] where p is a polynomial
    def __getitem__(self, m):
        return self.coef[-m - 1]

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
            c = [self.field.zero() for i in range(d1 - d2)] + other.coef
            return Polynomial(add_arrays(self.coef, c), self.field)
        else: 
            c = [self.field.zero() for i in range(d2 - d1)] + self.coef
            return Polynomial(add_arrays(other.coef, c), self.field)

    def __neg__(self):
        return Polynomial([-c for c in self.coef], self.field)

    def __sub__(self, other):
        assert self.field == other.field
        return self + (-other)

    # multiply this polynomial by (x ^ power)
    def multiply_monic_mononomial(self, degree):
        return Polynomial(self.coef + [self.field.zero() for i in range(degree)], self.field)

    def __mul__(self, other):
        if self.is_zero():
            return self
        elif type(other) is type(self.coef[0]):
            # in other words: multiply this polynomial by a scalar that is in the same field as the coefficients
            return self.scale(other)
        else:
            assert type(other) is Polynomial and other.field == self.field
            p, q = self, other
            # p(x) * (ax^2 + bx + c) = ax^2 * p(x) + bx * p(x) + c*p(x)
            prod = Polynomial.zero(self.field)
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
def divide_polynomials(first, second):
    assert first.field == second.field
    F = first.field
    d1 = first.degree()
    d2 = second.degree()
    if d1 < d2:
        return Polynomial.zero(F), first
    deg_quotient = d1 - d2

    # long division
    leading_coef = first.coef[0] / second.coef[0]
    quotient_leading_terms = leading_coef * Polynomial.monic_mononomial(deg_quotient, F)
    second_times_quot = (leading_coef * second).multiply_monic_mononomial(deg_quotient)
    smaller_polynomial = first - second_times_quot # The coefficient on the leading term will disappear
    if smaller_polynomial.is_zero():
        return quotient_leading_terms, Polynomial.zero(F)
    quotient_smaller_polynomial, remainder = divide_polynomials(smaller_polynomial, second)
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