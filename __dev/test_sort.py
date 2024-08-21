from itertools import product

import numpy as np

# # Given length n
# n = 3

# # Generate all binary tuples of length n
# binary_tuples = list(product([0, 1], repeat=n))

# # Convert tuples to a numpy array for easy manipulation
# binary_array = np.array(binary_tuples)
binary_array = np.array(
    [
        [0, 0, 1],
        [1, 0, 1],
        [0, 0, 1],
        [1, 0, 1],
        [0, 0, 1],
        [0, 0, 1],
        [1, 0, 1],
        [0, 0, 1],
    ],
    dtype=np.bool_,
)

# Use numpy's argsort to get the indices that would sort the array lexicographically

# this would be soooo much easier in C...
# have to try it later on.


@numba.njit
def to_int(array):
    res = np.uint32(0)
    for i, el in enumerate(array):
        res *= 2
        res += el
    return res


@numba.njit
def as_ints(matrix):
    results = np.empty(shape=(matrix.shape[0]), dtype=np.uint32)
    for i, row in enumerate(matrix):
        results[i] = to_int(row)
    return results

%%timeit
np.argsort(as_ints(binary_array))




# we could also count while doing all of that.

%%timeit
sort_order = np.argsort(
    binary_array.view([("", binary_array.dtype)] * binary_array.shape[1]), axis=0
).ravel()

np.argsort(binary_array, axis=0)

binary_array.view([("", binary_array.dtype)] * binary_array.shape[1])

# Display the original binary tuples and the sort order
print("Original binary tuples:")
print(binary_array)

print("\nSort order vector (numpy argsort equivalent):")
print(sort_order)

# Sorted binary tuples using the sort order
sorted_binary_tuples = binary_array[sort_order]

print("\nSorted binary tuples:")
print(sorted_binary_tuples)
# THIS WORKS: the trick is to use structured arrays that allows argsort to treat tuples as entries for sorting and likley sorting is implemented for them lexicographically by default.

a = np.array([[[1,0,1],
               [0,1,0]],
              [[1,1,0],
               [0,0,1]]])
b = np.packbits(a, axis=-1)
