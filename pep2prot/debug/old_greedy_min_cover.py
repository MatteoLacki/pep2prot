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
        # optimization: reporting supported nodes
        G.remove_nodes_from([n for n in covering() if G.degree(n) == 0])# non-covering after all
        # simplifying the problem by reporting nodes which already are uniquely supported
        cover = {n for m in to_cover() if G.degree(m)==1 for n in G[m]}
        # and removing all points next to them
        G.remove_nodes_from([m for n in cover for m in G[n]])
        G.remove_nodes_from(cover)

        # proper algorithm
        to_cover_cnt = G.B_cnt() if A_covers_B else G.A_cnt()
        while to_cover_cnt > 0:# this migh be a wrong criterion.
            n_max = max(covering(), key=lambda n: len(G[n]))# nodes with max degree
            print(n_max)
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