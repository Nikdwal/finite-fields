from util import add_arrays

# A field whose size is a prime number. This is isomorphic to the integers mod p.
class PrimeField:

    def __init__(self, p):
        assert p > 2 and not [i for i in range(2,p) if p % i == 0] # prime

        self.__size = p

        # define elements
        elems = []
        for i in range(p):
            elems.append(FieldElement(str(i), self))

        self._zero_elem = elems[0]
        self._one_elem = elems[1]

        # define addition
        self.add_table = {}
        for i in range(p):
            x = elems[i]
            for j in range(p):
                y = elems[j]
                self.add_table[x, y] = elems[(i + j) % p]

        # define symmetric elements
        self.negatives = {}
        for i in range(p):
            self.negatives[elems[i]] = elems[(p - i) % p]

        # find a multiplicative generator and its powers
        self._generator_powers = [None for i in range(p - 1)]    # generator_powers[i] = g^i
        self._generator_exponents = {} # generator_exponents[g^i] = i
        for i in range(1, p):
            generated_values = self._cycle_exp(i)
            if len(generated_values) == p - 1:
                # this is a generator of the multiplicative group
                for exp in range(p - 1):
                    self._generator_powers[exp] = elems[generated_values[exp]]
                    self._generator_exponents[elems[generated_values[exp]]] = exp
                break

        if not self._generator_exponents:
            raise ValueError("Error: expected to find a generator in GF(" + p + ").")

    # get the n'th element in the field
    # WARNING: the elements in a finite field do not follow a particular total order
    # the order used here is given by successive powers of the generator
    # use this as self[i] or F[i] when F == self
    def __getitem__(self, i):
        return self.get_elems()[i]

    # use this in for-loops
    # say self == F, then you can write: "for elem in F:..."
    def __iter__(self):
        return self.get_elems().__iter__()

    def __str__(self):
        return [str()]

    def __len__(self):
        return self.__size

    def zero(self):
        return self._zero_elem

    def one(self):
        return self._one_elem

    def generator(self):
        return self._generator_powers[1]

    # if the generator of the multiplicative group is elem = g^i, this returns i
    def generator_exponent(self, elem):
        return self._generator_exponents[elem]

    # returns g^exp where g is the generator of the multiplicative group
    def generator_to_power(self, exp):
        # the multiplicative group has size len(self) - 1
        return self._generator_powers[exp % (len(self) - 1)]

    def negative(self, elem):
        return self.negatives[elem]

    # find all the powers of a given (integer) value mod p
    def _cycle_exp(self, value):
        assert(type(value) is int)
        v = 1
        values = [v]
        for i in range(1, len(self)):
            v = (v * value) % len(self)
            if v in values:
                return values
            else:
                values.append(v)

    def get_elems(self):
        return self._generator_powers + [self.zero()]

class FieldElement:
    def __init__(self, name, field):
        self.name = name
        self.field = field
        self.__attrs = (self.name, self.field)

    def is_zero(self):
        return self == self.field.zero()

    def is_one(self):
        return self == self.field.one()

    def __eq__(self, other):
        return self.__attrs == other.__attrs

    def __hash__(self):
        return hash(self.__attrs)

    def __str__(self):
        return self.name

    def __add__(self, other):
        assert self.field == other.field
        return self.field.add_table[self, other]

    def __neg__(self):
        return self.field.negative(self)

    def __sub__(self, other):
        assert self.field == other.field
        return self + (-other)

    def __pow__(self, power, modulo=None):
        exp = self.field.generator_exponent(self)
        return self.field.generator_to_power(exp * power)

    def __mul__(self, other):

        if type(other) is FieldElement:
            assert self.field == other.field
            if self.field.zero() in [self, other]:
                return self.field.zero()

            e1 = self.field.generator_exponent(self)
            e2 = self.field.generator_exponent(other)
            return self.field.generator_to_power(e1 + e2)
        else:
            return other.scale(self)

    def __rmul__(self, n):
        # scale this element by n (add it to itself n times)
        assert type(n) is int
        s = self.field.zero()
        for i in range(n % len(self.field)):
            s += self
        return s

    # Divide self by other. You can always divide by 'other' because this is a field.
    def __truediv__(self, other):
        assert self.field == other.field
        exp = self.field.generator_exponent(other)
        inv = self.field.generator_to_power(-exp)
        return self * inv


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
        return Polynomial.x_to_power(0, field)

    @staticmethod
    def x_to_power(degree, field):
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
    def multiply_x_to_power(self, degree):
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
                prod += q[i] * p.multiply_x_to_power(i)
            return prod

    # multiply by a scalar
    def scale(self, scalar):
        assert type(scalar) is type(self.coef[0])
        coefs = [scalar * c for c in self.coef]
        return Polynomial(coefs, self.field)

    # return the quotient and the remainder
    def __divmod__(self, other):
        assert self.field == other.field
        return Polynomial.divide_polynomials(self, other)

    # return the quotient, but not the remainder
    def __floordiv__(self, other):
        assert self.field == other.field
        return divmod(self, other)[0]

    # return the remainder for self / other
    def __mod__(self, other):
        assert self.field == other.field
        return divmod(self, other)[1]

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
    def all_monic_polynomials(degree, field):
        if degree == 0:
            yield Polynomial([field.one()], field)
        else:
            for c in field:
                for p in Polynomial.all_monic_polynomials(degree - 1, field):
                    yield p.multiply_x_to_power(1) + Polynomial([c], field)

    @staticmethod
    def _make_larger_irreducible_polyns(field, irreducibles):
        # there are a few relatively clever tricks to find all of these given a certain field
        # however, none of them are quite remarkable; nor are they easy to implement algorithmically
        for p in irreducibles:
            for c in field:
                # this will make all monic polynomials of degree one higher
                w = p.multiply_x_to_power(1) + Polynomial([c], field)

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

    # if this is a primitive polynomial, return the successive powers of the primitive element
    # otherwise, return None
    def get_generated_polynomials(self):
        if self.is_reducible():
            return None

        F = self.field
        w = self
        q = len(F)
        k = w.degree()
        alpha = Polynomial.x_to_power(1, F)
        alpha_to_power_i = alpha
        generated_elems = [Polynomial.one(F), alpha]
        for i in range(2, q**k - 1):
            alpha_to_power_i = (alpha * alpha_to_power_i) % w
            if alpha_to_power_i in generated_elems:
                return None
            generated_elems.append(alpha_to_power_i)

        return generated_elems

    @staticmethod
    def find_primitive_polynomials(degree, field):
        for p in Polynomial.all_monic_polynomials(degree, field):
            polyns =  p.get_generated_polynomials()
            if polyns is not None:
                yield p, polyns

    @staticmethod
    def find_primitive_polynomial(degree, field):
        for p, polyns in Polynomial.find_primitive_polynomials(degree, field):
            return p, polyns

    # returns the quotient and the remainder
    @staticmethod
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
        quotient_leading_terms = leading_coef * Polynomial.x_to_power(deg_quotient, F)
        denom_times_quot = (leading_coef * denominator).multiply_x_to_power(deg_quotient)
        smaller_polynomial = numerator - denom_times_quot # The coefficient on the leading term will disappear
        if smaller_polynomial.is_zero():
            return quotient_leading_terms, Polynomial.zero(F)
        quotient_smaller_polynomial, remainder = divide_polynomials(smaller_polynomial, denominator)
        return quotient_leading_terms + quotient_smaller_polynomial, remainder


#===== testing
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
# print(p % ((p * q) + Polynomial.x_to_power(1, Z5)))
# print(p - q)
# print(p(four))
#
# factors = p.factor()
# for factor in factors:
#     print(factor)
#
# pol = Polynomial.x_to_power(0, Z5)
# for factor in factors:
#     pol *= factor
#
# print(pol)
# pol = Polynomial.find_primitive_polynomial(2, Z5)
# print(pol[0])