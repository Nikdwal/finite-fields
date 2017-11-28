class PrimeField:

    # A field whose size is a prime number. This is isomorphic to the integers mod p.
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


#===== testing
Z5 = PrimeField(5)
elems= Z5.get_elems()
g = Z5.generator()
zero = Z5.zero()
one = Z5.one()
two = one + one
three = one + two
four = one + three