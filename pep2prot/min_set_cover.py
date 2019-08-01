class SetNonCoverableException(Exception):
    pass

def greedy_minimal_cover(G, A_covers_B=True, copy_G=True):
    """Get the greedy approximation to the minimal set cover.
    
    Args:
        G (BiGraph): The Bigraph of interest (best a connected component thereof).
        A_covers_B (boolean): should we cover set B with elements from A, or the opposite.
        copy_G (boolean): should we work on a copy of G?
    Returns:
        list: A nodes that cover the B nodes
    """
    if copy_G:
        G = G.copy()
    to_cover, covering = (G.B, G.A) if A_covers_B else (G.A, G.B)
    non_coverable = [m for m in to_cover() if G.degree(m) == 0]
    if len(non_coverable) > 0:
        raise SetNonCoverableException("This set is impossible to cover: {}".format(str(non_coverable)))

    covering_cnt = G.A_cnt() if A_covers_B else G.B_cnt()
    if covering_cnt > 1:
        # optimization: reporting supproted nodes
        G.remove_nodes_from([n for n in covering() if G.degree(n) == 0])# non-covering after all
        # simplifying the problem by reporting nodes which already are uniquely supported
        cover = {n for m in to_cover() if G.degree(m)==1 for n in G[m]}
        # and removing all points next to them
        G.remove_nodes_from([m for n in cover for m in G[n]])
        G.remove_nodes_from(cover)

        # proper algorithm
        to_cover_cnt = G.B_cnt() if A_covers_B else G.A_cnt()
        while to_cover_cnt > 0:
            n_max = max(covering(), key=lambda n: len(G[n]))
            cover.add(n_max)
            n_max_neighbors = list(G[n_max])
            G.remove_nodes_from(n_max_neighbors)
            G.remove_node(n_max)
            G.remove_nodes_from([n for n in covering() if G.degree(n) == 0])# covers what n does
            to_cover_cnt -= len(n_max_neighbors)
        return cover
    else: 
        return {next(covering())}


def test_greedy_minimal_cover():
    from pep2prot.graphs import BiGraph
    G = BiGraph([('a',0),('a',1),('b',1),('b',2),('c',2),('c',3),
                 ('d',3),('d',4),('e',4),('e',5),('f',6),('f',7),
                 ('g',7),('g',8),('g',9),('h',9),('h',10),('i',10),
                 ('i',11),('j',12),('j',14),('k',12),('k',13),
                 ('l',13),('l',14)])
    set.union(*(greedy_minimal_cover(cc) for cc in G.components()))

    # TEST   0-a-1-b-2-c-3 >> {a,c} cover it.  
    G = BiGraph([('a',0),('a',1),('b',1),('b',2),('c',2),('c',3)])
    assert greedy_minimal_cover(G) == {'a', 'c'}


#TODO: something is wrong here.
def greedy_minimal_cover_2(G, A_covers_B=True):
    """Get the greedy approximation to the minimal set cover.
    
    Args:
        G (BiGraph): The Bigraph of interest (best a connected component thereof).
        A_covers_B (boolean): should we cover set B with elements from A, or the opposite.
        copy_G (boolean): should we work on a copy of G?
    Returns:
        list: A nodes that cover the B nodes
    """
    G = G.copy()
    to_cover, covering = (G.B, G.A) if A_covers_B else (G.A, G.B)

    non_coverable = [m for m in to_cover() if G.degree(m) == 0]
    if len(non_coverable) > 0:
        raise SetNonCoverableException("This set is impossible to cover: {}".format(str(non_coverable)))

    covering_cnt = G.A_cnt() if A_covers_B else G.B_cnt()
    if covering_cnt > 1:
        cover = {n for m in to_cover() if G.degree(m)==1 for n in G[m]}# supported nodes
        G.remove_nodes_from([m for n in cover for m in G[n]])# remove what they cover
        G.remove_nodes_from(cover)# remove supported nodes
        G.remove_nodes_from([n for n in covering() if G.degree(n) == 0])# remove unsupported
        for C in G.components():
            cover.update(inner_greedy(C, A_covers_B))# recursive move
        return cover
    else: 
        return {next(covering())}


def inner_greedy(C, A_covers_B=True):
    """
    Args:
        C (BiGraph): Connected BiGraph.
        A_covers_B (boolean): should we cover set B with elements from A, or the opposite.
    Returns:
        set: covering nodes.
    """
    C = C.copy() # is this necessary, if it is inside????
    to_cover, covering = (C.B, C.A) if A_covers_B else (C.A, C.B)
    max_covering_edge_cnt = max(len(C[n]) for n in covering())
    cover = {n for n in covering() if len(C[n]) == max_covering_edge_cnt}
    covered = {m for n in cover for m in C[n]}
    C.remove_nodes_from(covered)
    C.remove_nodes_from(cover)
    C.remove_nodes_from([n for n in covering() if C.degree(n) == 0])# covering what is already covered
    for cc in C.components():
        cover.update(inner_greedy(cc, A_covers_B=A_covers_B))
    return cover


def test_greedy_minimal_cover_2():
    from pep2prot.graphs import BiGraph
    T = BiGraph([('A',1),('A',2),('E',1),('B',2),('B',3),('F',3),('B',4),
                ('C',4),('C',5),('G',5),('C',6),('I',6),('C',7),
                ('D',7),('D',8),('H',8),('D',9),('J',9),('D',10),('K',10)])
    # T.draw(with_labels=True)
    assert greedy_minimal_cover_2(T) == {'A','B','C','D'}

    T = BiGraph([('a',0),('a',1),('b',1),('b',2),('c',2),('c',3)])
    assert greedy_minimal_cover_2(T) == {'a','c'}

    T = BiGraph([('a',0),('a',1),('b',1),('b',2),('c',2),('c',3),
                 ('d',3),('d',4)])
    assert greedy_minimal_cover_2(T) == {'a','b','c','d'}

    # Interesing case:
    TSC = BiGraph([('A',1),('A',2),('A',3),
                   ('B',5),('B',2),('B',3),('B',4),
                   ('C',4),('C',5),
                   ('D',6),('D',7),('D',8),('D',9),
                   ('E',8),('E',9),
                   ('F',10),('F',12),
                   ('G',11),('G',12),
                   ('H',10),('H',11),('H',12)])
    TSC = TSC.form_groups()
    TSC_MC = TSC.greedy_minimal_cover()
    assert TSC_MC=={frozenset({'H'}), frozenset({'D'}), frozenset({'A'}), frozenset({'C'}), frozenset({'B'})}
    node_sizes = [40 if n in TSC_MC else 10 for n in TSC]
    # TSC.draw(node_size=node_sizes)



