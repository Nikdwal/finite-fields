class Polynomial:
    def __init__(self, coef, q):
        first_nonzero = -1
        coef = [round_elem(x) % q for x in coef]
        for i in range(len(coef)):
            if coef[i] != 0:
                first_nonzero = i
                break
        # print(coef, first_nonzero)
        self.coef = coef[first_nonzero:] if first_nonzero >= 0 else []
        self.q = q
        # if coef == np.array([1, 2, 1]) this represents x^2 + 2x + 1

    def __str__(self):
        if len(self.coef) == 0: 
            return "0"
        s = ""
        d = self.degree()
        sup = str.maketrans("0123456789", "⁰ ²³⁴⁵⁶⁷⁸⁹")
        for i in range(d):
            s += (str(self.coef[i]) if self.coef[i] != 1 else "") + "x" + str(d - i).translate(sup) + " + " 
        s += str(self.coef[d])
        return s

    def is_zero(self):
        return self.coef == [] 

    def degree(self):
        return len(self.coef) - 1
    
    def scale(self, scalar):
        return Polynomial(scale(self.coef, scalar), self.q)

    def add(self, other):
        assert self.q == other.q
        d1 = self.degree()
        d2 = other.degree()
        if d1 > d2:
            c = [0 for i in range(d1 - d2)] + other.coef
            return Polynomial(add_arrays(self.coef, c), self.q)
        else: 
            c = [0 for i in range(d2 - d1)] + self.coef
            return Polynomial(add_arrays(other.coef, c), self.q)

    def subtract(self, other):
        return self.add(other.scale(-1))

#    def multiply(self, other):
#        assert self.q == other.q
#        return Polynomial(np.convolve(self.coef, other.coef), self.q)

def round_elem(x):
    # todo: revise when working algebraically
    return int(x)

def scale(lst, scalar):
    return [scalar * x for x in lst]

def add_arrays(lst1, lst2):
    assert len(lst1) == len(lst2)
    s = [x for x in lst1]
    for i in range(len(lst1)):
        s[i] += lst2[i]
    return s

def divide(first, second):
    #print(first)
    #print("divided by ", second)
    # staartdeling
    assert first.q == second.q
    q = first.q
    d1 = first.degree()
    d2 = second.degree()
    if d1 < d2:
        return None
    deg_quotient = d1 - d2 
    leading_coef = divmodq(first.coef[0], second.coef[0], q)
    if leading_coef is None:
        return None
    p = Polynomial([leading_coef] + [0 for i in range(deg_quotient)], q)
    a = scale(second.coef, leading_coef) + [0 for i in range(deg_quotient)]  
    b = Polynomial(a, q)
    c = first.subtract(b)
    if c.is_zero():
        return p
    rec = divide(c, second)
    if rec == None:
        return None
    else:
        return p.add(rec)

def divmodq(a, b, q):
    a = a % q
    for i in range(1, q):
        if (i*b) % q == a:
            return i
    return None

u = Polynomial([1, 1, 2], 3)
v = Polynomial([2, 1], 3)
import numpy as np
w = Polynomial(np.convolve(u.coef, v.coef), 3)
divide(w, v)

