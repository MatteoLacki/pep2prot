def asserting_unique(x):
    x = set(x)
    assert len(x) == 1
    return x.pop()
