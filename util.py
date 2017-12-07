
def add_arrays(lst1, lst2):
    assert len(lst1) == len(lst2)
    s = [x for x in lst1]
    for i in range(len(lst1)):
        s[i] += lst2[i]
    return s

def product(l):
    # cannot start from prod = 0 because we need the zero element from a finite field
    prod = l[0]
    for i in range(1, len(l)):
        prod *= l[i]
    return prod

def dot_product(arr1, arr2):
    assert len(arr1) == len(arr2)
    # cannot start from dp = 0 because we need the zero element from a finite field
    dp = arr1[0] * arr2[0]
    for i in range(1, len(arr1)):
        dp += arr1[i] * arr2[i]
    return dp

greek_alphabet = "αβγδεζηθικλμνξοπρςστυφχψω"