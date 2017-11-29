
def add_arrays(lst1, lst2):
    assert len(lst1) == len(lst2)
    s = [x for x in lst1]
    for i in range(len(lst1)):
        s[i] += lst2[i]
    return s

def product(l):
    prod = l[0]
    for i in range(1, len(l)):
        prod *= l[i]
    return prod

greek_alphabet = "αβγδεζηθικλμνξοπρςστυφχψω"