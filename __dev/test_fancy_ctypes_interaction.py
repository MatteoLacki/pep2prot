import ctypes

import numpy as np

# Step 1: Create a NumPy array of uint32 elements
N = 10  # Size of the array
arr = np.zeros(N, dtype=np.uint32)

# Step 2: Expose the NumPy array's data as a ctypes pointer
# Get the pointer to the underlying data of the NumPy array
arr_ptr = arr.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))


# Step 3: Manipulate individual bits using the ctypes pointer
def set_bit(index, bit_position):
    """Set a specific bit to 1 in the uint32 element at index."""
    arr_ptr[index] |= 1 << bit_position


def clear_bit(index, bit_position):
    """Clear a specific bit to 0 in the uint32 element at index."""
    arr_ptr[index] &= ~(1 << bit_position)


def toggle_bit(index, bit_position):
    """Toggle a specific bit in the uint32 element at index."""
    arr_ptr[index] ^= 1 << bit_position


def check_bit(index, bit_position):
    """Check if a specific bit is set in the uint32 element at index."""
    return (arr_ptr[index] >> bit_position) & 1


# Example Usage
set_bit(0, 5)  # Set the 5th bit of the first element
toggle_bit(0, 7)  # Toggle the 7th bit of the first element
clear_bit(0, 5)  # Clear the 5th bit of the first element

# Check the bit value
print("Bit 7 of arr[0] is:", check_bit(0, 7))
print("Bit 5 of arr[0] is:", check_bit(0, 5))

# Print the modified array
print("Array:", arr)

X = np.array([[0, 1], [1, 0]], dtype=np.bool_)
X.view(dtype=np.int8)


x = np.array([(1, 2)], dtype=[("a", np.int8), ("b", np.int8)])
x.view(dtype=np.int8)
x.view(dtype=np.bool_)
# Viewing array data using a different type and dtype:
y = x.view(dtype=np.int16)
# this is very interesting: pure casting.
y
x

# Creating a view on a structured array so it can be used in calculations
x = np.array([(1, 2), (3, 4)], dtype=[("a", np.int8), ("b", np.int8)])
xv = x.view(dtype=np.int8).reshape(-1, 2)

# so we get a getitem interface for free now
xv[0, 1] = 20
x
