import numba
import numpy as np
import numpy.typing as npt


@numba.njit
def check_arrays_equal(a, b):
    assert len(a) == len(b)
    for i in range(len(a)):
        if a[i] != b[i]:
            return False
    return True


def get_dtype_and_its_bits():
    for x in (32, 64):
        if str(x) in str(np.uintp):
            return np.uintp, x


@numba.njit
def check_vecs_the_same(array: npt.NDArray) -> bool:
    if len(array) <= 1:
        return True
    prev = array[0]
    for i in range(len(array)):
        if not np.all(prev == array[i]):
            return False
        prev = array[i]
    return True
