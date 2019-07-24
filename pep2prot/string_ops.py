"""Here, we establish all indices which a substring p occupies in a string psp.

There are two possible approaches: either we count only disjoin occurences,
so searching for 'aaa' in 'aaaaaa' would result in (0,3),(3,6), or all occurences,
resulting in (0,3),(1,4),(2,5),(3,6).

For peptides it makes sense to look for the first case, as it is unlikely to meet 
situations like 'aaa'.
"""

import re


def starts_and_ends(p, psp):
    for m in re.finditer("(?=" + p + ")", psp):
        s = m.start()
        e = s + len(p)
        yield (s,e)

def starts_and_ends2(p, psp):
    for m in re.finditer(p, psp):
        s = m.start()
        e = s + len(p)
        yield (s,e)

def findall(p, psp):
    w = len(p)
    i = psp.find(p)
    while i != -1:
        yield (i,i+w)
        i = psp.find(p, i+1)

def find_all(p, psp):
    start = 0
    w = len(p)
    while True:
        start = psp.find(p, start)
        if start == -1: return
        yield (start,start+w)
        start += len(p) # use start += 1 to find overlapping matches

def find_all2(p, psp):
    start = 0
    w = len(p)
    while True:
        start = psp.find(p, start)
        if start == -1: return
        yield (start, start+w)
        start += 1

def find_all3(p, psp):
    start = 0
    w = len(p)
    res = []
    while True:
        start = psp.find(p, start)
        if start == -1: return res
        res.append((start, start+w))
        start += 1

def find_all4(p, psp):
    start = 0
    w = len(p)
    res = []
    while True:
        start = psp.find(p, start)
        if start == -1: return res
        res.append((start, start+w))
        start += len(p)

def find_indices(p, psp):
    w = len(p)
    i = 0
    x = psp.split(p)
    i += len(x[0])
    yield i, i+w
    i += w
    for h in x[1:-1]:
        i += len(h)
        yield i, i+w
        i += w

def find_indices2(p, psp):
    w = len(p)
    i = 0
    for h in psp.split(p)[0:-1]:
        i += len(h)
        yield i, i+w
        i += w

def find_indices3(p, psp):
    if p:
        w = len(p)
        i = 0
        res = []
        for h in psp.split(p)[0:-1]:
            i += len(h)
            res.append((i,i+w))
            i += w
        return res
    else:
        return []

# p = 'MQIFVK'
# psp = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLDVEPSVTTKKVKQEDRRTFLTTVSKKSPPCACSWV'
# %%timeit
# list(starts_and_ends(p,psp))
# %%timeit
# list(starts_and_ends2(p,psp))
# %%timeit
# list(findall(p,psp))
# %%timeit
# list(find_all(p,psp))
# %%timeit
# list(find_all2(p,psp))
# %%timeit
# find_all3(p,psp)
# %%timeit
# find_all4(p,psp)
# %%timeit
# list(find_indices(p,psp))
# %%timeit
# list(find_indices2(p,psp))
# %%timeit
# find_indices3(p,psp)
