from util import *
from abc import ABC
from math import gcd

class FiniteField(ABC):
    '''
    Any instance of FiniteField must at least have the following variables:
        - _elems_by_value: a dictionary that gives a fieldElem for each value
        - _zero_elem
        - _one_elem
        - _generator_powers: if the generator is g, then _generator_powers[i] must be g[i]
        - _generator_exponents: a dict that gives i for every g[i]
        - divisor: if this group is S mod m, then m == divisor, for example: the divisor of Z mod 5 is 5, the divisor of GF(q) mod w is w
    And call the following methods after construction:
        - _make_elems_by_repr
    '''

    def __eq__(self, other):
        # TODO: This should be enough, right? In principle, all fields with the same size are isomorphic.
        # however, this does not mean that they are 'equal'
        # for instance, this method differentiates between two different extended fields with a different divisor polynomial
        return isinstance(other, FiniteField) \
               and len(self) == len(other) \
               and self.divisor == other.divisor \
               and self.generator().name == other.generator().name

    def __hash__(self):
        return hash((len(self), self.divisor))

    def is_subfield(self, other):
        return self <= other

    # strict subfield
    def __lt__(self, other):
        return isinstance(other, ExtendedField) and self <= other.subfield

    # non-strict subfield
    def __le__(self, other):
        return self == other or self < other

    def  __gt__(self, other):
        return other < self

    def __ge__(self, other):
        return other <= self

    # return the "largest" of self and other
    # largest means that the one that is returned is a superfield of the other one
    def largest(self, other):
        if self <= other:
            return other
        elif other <= self:
            return self
        else:
            return None

    # make a dictionary of all the elements in this field, with their string (not repr) representations as keys
    def _make_elems_by_repr(self):
        self._elems_by_repr = {}
        for elem in self:
            self._elems_by_repr[str(elem)] = elem

    # given the representations of one or more elements (e.g. a number or a string), return the elements themselves
    # This can return a KeyError if the element is not found
    # reprs could be a singleton, a list, or a list of lists
    # eg use as F[1,2,"α"]
    def __getitem__(self, *items):
        reprs = items[0]

        if isinstance(reprs, int):
            # go back to the elementary field that this is field is based on
            # (an IntegerField) and get the item for reprs % p
            # e.g. GF_5[7] should return the FieldElement 2.
            base_field = self
            while isinstance(base_field, ExtendedField):
                base_field = base_field.subfield
            # base_field should now be an IntegerField
            return self[base_field.get_elem_by_value(reprs)]

        if isinstance(reprs, list):
            # This is a list. Return a list.
            return [self[bag] for bag in reprs]
        elif isinstance(reprs, tuple):
            # This is a tuple. Return a tuple
            return tuple([self[bag] for bag in reprs])
        else:
            # This is a singleton. Return a singleton
            return self._elems_by_repr[str(reprs)]

    # F.cast(x) should always be equal to x when x is an element of F.
    def cast(self, item):
        return self[item]

    def __call__(self, *args, **kwargs):
        return self.__getitem__(args)

    # use this in for-loops
    # say self == F, then you can write: "for elem in F:..."
    def __iter__(self):
        return self.get_elems().__iter__()

    def __repr__(self):
        return str(self.get_elems())

    def __len__(self):
        return len(self._elems_by_value)

    def zero(self):
        return self._zero_elem

    def one(self):
        return self._one_elem

    def generator(self):
        if len(self) == 2:
            # the generator is one
            return self.one()
        return self._generator_powers[1]

    # if the generator of the multiplicative group is elem = g^i, this returns i
    def generator_exponent(self, elem):
        if elem.is_zero():
            raise ZeroDivisionError("Cannot find the exponent for zero. The generator is not a zero divisor.")
        return self._generator_exponents[elem]

    # returns g^exp where g is the generator of the multiplicative group
    def generator_to_power(self, exp):
        # the multiplicative group has size len(self) - 1
        return self._generator_powers[exp % (len(self) - 1)]

    def get_elems(self):
        # the method "values" is misleading, because it comes from the dict class. It actually refers to the elements themselves.
        return list(self._elems_by_value.values())

    def get_elem_by_value(self, val):
        return self._elems_by_value[val % self.divisor]


    # ==== start polynomial-related methods
    def x(self):
        return Polynomial.x_to_power(1, self)

    def all_monic_polynomials(self, degree):
        if degree == 0:
            yield Polynomial([self.one()])
        else:
            for c in self:
                for p in self.all_monic_polynomials(degree - 1):
                    yield p.multiply_x_to_power(1) + Polynomial([c])

    def _make_larger_irreducible_polyns(self, irreducibles):
        # there are a few relatively clever tricks to find all of these given a certain field
        # however, none of them are quite remarkable; nor are they easy to implement algorithmically
        for p in irreducibles:
            for c in self:
                # this will make all monic polynomials of degree one higher
                w = p.multiply_x_to_power(1) + Polynomial([c])

                if not w.is_reducible():
                    yield w

    def find_primitive_polynomials(self, degree):
        for p in self.all_monic_polynomials(degree):
            polyns = p._get_generated_polynomials()
            if polyns is not None:
                yield p

    def find_primitive_polynomial(self, degree):
        for p in self.find_primitive_polynomials(degree):
            return p

    def nth_root_of_unity(self, n):
        q = len(self)
        if gcd(q, n) == 1:
            for k in range(1, n+1):
                if (q**k - 1) % n == 0:
                    l = (q**k - 1) // n
                    extField = ExtendedField(self, k, "α")
                    return extField.generator_to_power(l)
        else:
            # the usual method won't work here
            raise NotImplementedError()

    # factor x^n - 1
    def factor_nth_root(self, n):
        q = len(self)
        GFq = self
        beta = self.nth_root_of_unity(n)
        cyclo_cosets = cyclotomic_cosets(q, n)
        factors_in_GFqk = []
        GFqk = beta.field
        for coset in cyclo_cosets:
            factor = Polynomial.one(GFqk)
            for exp in coset:
                factor *= Polynomial([-beta**exp, GFqk.one()])
            factors_in_GFqk.append(factor)
        # Each of the factors_in_GFqk is now guaranteed to be a polynomial in GF(q)[x]
        # but their data types are polynomials in GF(q^k)[x] (hence the name of the variable)
        # This effectively means that each factor is Sum(p_i(alpha) * x^i) where p_i is a constant polynomial in GF(q).
        # Which means their data type is a Polynomial whose coefficients are FieldElements whose
        # "value" attribute is this constant p_i(alpha), from which we need to extract the coefficient on the constant term
        factors_in_GFq = []
        for factor in factors_in_GFqk:
            factors_in_GFq.append(Polynomial([coefficient.value[0] for coefficient in factor]))
        return factors_in_GFq

    # === end polynomial-related methods

class IntegerField(FiniteField):

    def __init__(self, p):
        if not (p >= 2 and not [i for i in range(2,p) if p % i == 0]):
            raise ValueError("The size of Z mod p must be prime. Otherwise it would not be a field.")

        self.divisor = p

        # define elements
        elems = []
        for i in range(p):
            elems.append(FieldElement(i, str(i), self))

        self._zero_elem = elems[0]
        self._one_elem = elems[1]

        # this is a dictionary because in general, fields are not restricted to numbers
        self._elems_by_value = {}
        for i in range(p):
            self._elems_by_value[i] = elems[i]

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

        # TODO: is finding a generator guaranteed (in other words: is this check redundant?)
        if not self._generator_exponents:
            raise ValueError("Error: expected to find a generator in GF(" + p + ").")

        self._make_elems_by_repr()

    # find all the powers of a given (integer) value mod p
    def _cycle_exp(self, value):
        assert isinstance(value, int)
        v = 1
        values = [v]
        for i in range(1, len(self)):
            v = (v * value) % len(self)
            if v in values:
                return values
            else:
                values.append(v)

# Makes an extended field, given a subfield
# The most important thing to remember about its internal workings is that each element in this field
# is formally a FieldElem whose "value" property is a polynomial of degree k over the smaller field
class ExtendedField(FiniteField):
    def __init__(self, subfield, k, name_primitive_elem, primitive_polynomial=None):
        if primitive_polynomial is None:
            w = subfield.find_primitive_polynomial(k)
        else:
            w = primitive_polynomial
        quotient_modules = w._get_generated_polynomials() # "restklassen mod w"

        self.divisor = w
        self._generator_exponents = {}
        self._elems_by_value = {}

        zero_polynomial = Polynomial.zero(subfield)
        for p in [zero_polynomial] + quotient_modules:
            self._elems_by_value[p] = FieldElement(p, str(p).replace("x", name_primitive_elem), self)

        self._generator_powers = [self._elems_by_value[val] for val in quotient_modules]
        for i in range(len(quotient_modules)):
            self._generator_exponents[self._generator_powers[i]] = i

        self._zero_elem = self._elems_by_value[zero_polynomial]
        self._one_elem = self._generator_powers[0]

        self._make_elems_by_repr()

        self.subfield = subfield

class FieldElement:
    def __init__(self, value, name, field):
        self.value = value
        self.name = name
        self.field = field
        self.__attrs = (self.value, self.name, self.field)

    def is_zero(self):
        return self == self.field.zero()

    def is_one(self):
        return self == self.field.one()

    def __eq__(self, other):
        return isinstance(other, FieldElement) and self.__attrs == other.__attrs

    def __hash__(self):
        return hash(self.__attrs)

    def __repr__(self):
        return self.name

    # add two elements with automatic type conversion
    def __add__(self, other):
        if isinstance(other, FieldElement):
            F = self.field.largest(other.field)
            if F is None:
                raise ValueError("Cannot add these two elements. Their fields are incompatible.")
            val1 = F.cast(self).value
            val2 = F.cast(other).value
            return F.get_elem_by_value(val1 + val2)

        elif isinstance(other, Polynomial):
            p = Polynomial(self)
            # the following will raise an exception if we cannot add p to other
            return p + other

        else: # usually an integer
            # this may raise a KeyError if other is not in the field
            val1 = self.value
            val2 = self.field.cast(other).value
            return self.field.get_elem_by_value(val1 + val2)

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        return self.field.get_elem_by_value(- self.value)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        # cast other to something we can subtract self from
        other = self.field.zero() + other
        return other - self

    def __pow__(self, power, modulo=None):
        if self.is_zero():
            return self
        exp = self.field.generator_exponent(self)
        return self.field.generator_to_power(exp * power)

    def __mul__(self, other):

        if isinstance(other, FieldElement):
            F = self.field.largest(other.field)
            if F is None:
                raise ValueError("Cannot multiply these two elements. Their fields are incompatible.")
            elem1 = F.cast(self)
            elem2 = F.cast(other)

            if self.is_zero() or other.is_zero():
                return F.zero()

            # compute product using generators
            e1 = F.generator_exponent(elem1)
            e2 = F.generator_exponent(elem2)
            return F.generator_to_power(e1 + e2)

        elif isinstance(other, int):
            return self * self.field.cast(other)
        else:
            return other.scale(self)

    def __rmul__(self, other):
        return self * other

    def invert(self):
        exp = self.field.generator_exponent(self)
        return self.field.generator_to_power(-exp)

    # Divide self by other. You can always divide by 'other' because this is a field.
    def __truediv__(self, other):
        if isinstance(other, FieldElement):
            return self * other.invert()
        elif isinstance(other, int):
            return self / self.field.cast(other)

    def __rtruediv__(self, other):
        return self.invert() * other

class Polynomial:

    # coef is a list of coefficients from lowest to highest degree terms
    # or it could be a single field element
    def __init__(self, coef):
        if isinstance(coef, FieldElement):
            coef = [coef]
        else:
            assert(isinstance(coef, list) or isinstance(coef, tuple))
        # all elements must be in the same field
        assert all([c.field == coef[0].field for c in coef])

        # infer the field from the coefficients
        self.field = coef[0].field

        # remove higher degree terms if their coefficients are zero (those would be redundant)
        first_zero = len(coef)
        for i in range(len(coef) -1, -1, -1):
            if coef[i].is_zero():
                first_zero -= 1
            else:
                break
        self.coef = list(coef[:first_zero])

        # If all elements are zero, self.coef is empty. However, this is rather difficult to work with.
        if not self.coef:
            self.coef = [self.field.zero()]

    @staticmethod
    def zero(field):
        return Polynomial([field.zero()])

    @staticmethod
    def one(field):
        return Polynomial.x_to_power(0, field)

    @staticmethod
    def x_to_power(degree, field):
        return Polynomial([field.zero() for i in range(degree)] + [field.one()])

    def __repr__(self):
        if self.is_zero():
            return str(self.field.zero())

        superscript = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
        s = "" if self[0].is_zero() else str(self[0])
        for i in range(1, len(self)):
            if not self.coef[i].is_zero():
                s += ("" if s == "" else " + ") +\
                     ( "" if self[i].is_one() else
                       ("(" + str(self[i]) + ")" if " " in str(self[i]) else str(self[i])) )+\
                     "x" + ("" if i == 1 else str(i).translate(superscript))
        return s

    # compute the value of the polynomial when substituting x=args[0]
    def __call__(self, *args, **kwargs):
        x = args[0]
        if isinstance(x, FieldElement):
            F = x.field.largest(self.field)
            if F is None:
                raise ValueError("Cannot apply this polynomial to this argument. The fields are incompatible")
            x = F.cast(x)
            p = self.cast_to_field(F)

            res = p[-1]
            for i in range(p.degree() - 1, -1, -1):
                res = p[i] + res * x
            return res
        else:
            # try casting it to a field element
            return self(self.field.cast(x))

    # get the n'th degree coefficient
    # usage: p[n] where p is a polynomial
    def __getitem__(self, m):
        if m >= len(self):
            # all higher order coefficients are zero, so just return that
            return self.field.zero()
        return self.coef[m]

    def __iter__(self):
        return self.coef.__iter__()

    def __eq__(self, other):
        return isinstance(other, Polynomial) and self.field == other.field and (self - other).is_zero()

    def __hash__(self):
        # TODO: is this correct? does a == b => hash(a) == hash(b) in finite fields?
        return hash(tuple(self.coef))

    def is_zero(self):
        return self.coef == [self.field.zero()]

    def is_one(self):
        return len(self) == 1 and self.coef[0].is_one()

    def is_power_of_x(self):
        return self[-1].is_one() and all([self[i].is_zero() for i in range(self.degree())])

    def __len__(self):
        return len(self.coef)

    def degree(self):
        return len(self) - 1

    def __add__(self, other):
        if isinstance(other, FieldElement):
            return other + self
        if isinstance(other, int):
            return self.field.cast(other) + self
        if isinstance(other, Polynomial):
            F = self.field.largest(other.field)
            if F is None:
                raise ValueError("Cannot add these two polynomials. Their fields are incompatible.")
            p1 = self.cast_to_field(F)
            p2 = other.cast_to_field(F)

            d1 = p1.degree()
            d2 = p2.degree()
            if d1 > d2:
                c = p2.coef + [F.zero() for i in range(d1 - d2)]
                return Polynomial(add_arrays(p1.coef, c))
            else:
                c = p1.coef + [F.zero() for i in range(d2 - d1)]
                return Polynomial(add_arrays(p2.coef, c))
        else:
            raise TypeError("Addition between Polynomial and " + str(type(other)) + " is not supported.")

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        return Polynomial([-c for c in self.coef])

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return (-self) + other

    # multiply this polynomial by (x ^ power)
    def multiply_x_to_power(self, degree):
        return Polynomial([self.field.zero() for i in range(degree)] + self.coef)

    def __mul__(self, other):
        if isinstance(other, FieldElement) or isinstance(other, int):
            return self.scale(other)
        elif isinstance(other, Polynomial):
            p, q = self, other
            # p(x) * (q_0 + q_1 x + q_2 x^2) = q_0 * p(x) + q_1 x * p(x) + q_2 x^2 * p(x)
            prod = Polynomial.zero(self.field)
            for i in range(len(other)):
                # this will promote automatically
                prod += q[i] * p.multiply_x_to_power(i)
            return prod
        else:
            raise TypeError("Multiplication between Polynomial and " + str(type(other)) + " is not supported.")

    def __rmul__(self, other):
        return self * other

    # multiply by a scalar
    def scale(self, scalar):
        if isinstance(scalar, FieldElement):
            F = scalar.field.largest(self.field)
            if F is None:
                raise ValueError("Cannot scale this polynomial by this scalar. The types are incompatible")
            scalar = F.cast(scalar)
            p = self.cast_to_field(F)
            return Polynomial([scalar * c for c in p.coef])
        elif isinstance(scalar, int):
            return self.scale(self.field.cast(scalar))
        else:
            raise TypeError("Scaling between" + str(type(scalar)) + " and Polynomial is not supported.")

    def __pow__(self, n, modulo=None):
        # first check if this is a power of x (as is usually the case when this method is called)
        if self.is_power_of_x():
            # this is more efficient than the general case
            m = self.degree()
            return self.multiply_x_to_power(m*n - m)
        else:
            res = Polynomial.one(self.field)
            for i in range(n):
                res *= self
            return res

    # return the quotient and the remainder
    def __divmod__(self, other):
        if isinstance(other, Polynomial):
            F = self.field.largest(other.field)
            if F is None:
                raise ValueError("Cannot divide these two polynomials. The remainder is non-zero.")
            p1 = self.cast_to_field(F)
            p2 = other.cast_to_field(F)
            return Polynomial._divide_polynomials(p1, p2)
        else:
            raise TypeError("Divmod between Polynomial and " + str(type(other)) + " is not supported.")

    # return the quotient, but not the remainder
    def __floordiv__(self, other):
        return divmod(self, other)[0]

    # return the remainder for self / other
    def __mod__(self, other):
        return divmod(self, other)[1]

    # return the quotient iff the remainder is zero
    def __truediv__(self, other):
        if isinstance(other, FieldElement) or isinstance(other, int):
            return self.scale(other ** -1)

        quotient, remainder = divmod(self, other)
        if remainder.is_zero():
            return quotient
        else:
            raise ValueError("Cannot divide these two polynomials. The remainder is non-zero.")

    def differentiate(self):
        new_coefs = self.coef[1:]
        for i in range(len(new_coefs)):
            new_coefs[i] = (i+1) * new_coefs[i]
        if len(new_coefs) == 0:
            # self is a constant polynomial:
            new_coefs = [self.field.zero()]
        return Polynomial(new_coefs)

    # find an x such that p(x) is zero
    def find_roots(self):
        roots = []
        # there is no clever way to do this in finite fields
        for x in self.field:
            if self(x).is_zero():
                roots.append(x)
        return roots

    def factor(self):
        # TODO: this part is very messy and inefficient
        F = self.field

        a_n = self[-1]
        if a_n.is_one():
            factors = []
            w = self
        else:
            # split off the leading coefficient
            factors = [Polynomial([a_n])]
            w = (F.one() / a_n) * self

        irreducibles = [Polynomial([F.one()])]
        for deg in range(1, 1 + self.degree() // 2):
            for p in F._make_larger_irreducible_polyns(irreducibles):

                if w.is_one():
                    # completely factorised
                    return factors

                quot, rem = divmod(w, p)
                if rem.is_zero():
                    # found a factor
                    factors.append(p)
                    w = quot

        # TODO: will this line ever be executed?
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

        irreducibles = [Polynomial([F.one()])]
        for deg in range(1, 1 + self.degree() // 2):
            for p in F._make_larger_irreducible_polyns(irreducibles):
                quot, rem = divmod(w, p)
                if rem.is_zero():
                    # found a factor
                    return True
        return False

    # if this is a primitive polynomial, return the successive powers of the primitive element
    # otherwise, return None
    def _get_generated_polynomials(self):
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

    # returns the quotient and the remainder
    @staticmethod
    def _divide_polynomials(numerator, denominator):
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
        quotient_smaller_polynomial, remainder = Polynomial._divide_polynomials(smaller_polynomial, denominator)
        return quotient_leading_terms + quotient_smaller_polynomial, remainder

    # self.cast_to_field(self.field) should always return self.
    def cast_to_field(self, field):
        # TODO: catch KeyError
        new_coefs = [field.cast(coef) for coef in self]
        return Polynomial(new_coefs)

# cyclotomic cosets mod n on GF(q)
def cyclotomic_cosets(q, n):
    cosets = []
    generated = [False for i in range(n)]
    for s in range(n):
        if generated[s]:
            continue
        coset = [s]
        sq_i = s
        while True:
            generated[sq_i] = True
            sq_i = (sq_i * q) % (n)
            if sq_i == s:
                break
            coset.append(sq_i)
        cosets.append(coset)
    return cosets


#===== testing
if __name__ == "__main__":
    pass