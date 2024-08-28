import typing


def asserting_unique(x):
    x = set(x)
    assert len(x) == 1
    return x.pop()


def is_nondecreasing(xx: typing.Iterable) -> bool:
    xx = iter(xx)
    x_prev = next(xx)
    for x in xx:
        if x < x_prev:
            return False
        x_prev = x
    return True
