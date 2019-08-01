def greedy_minimal_cover(G):
    G = G.copy()
    cover = set([])
    if G.prot_cnt() > 1:
        cover.update(r for p in G.peps() if G.degree(p)==1 for r in G[p])# supported proteins
        G.remove_nodes_from([p for r in cover for p in G[r]])# remove peptides they cover
        G.remove_nodes_from(cover)# remove supported proteins
        G.remove_nodes_from([r for r in G.prots() if G.degree(r) == 0])# remove unsupported prots
        if G.prot_cnt() > 1:
            max_cover_degree = max(G.degree(r) for r in G.prots())
            max_degree_nodes = {r for r in G.prots() if G.degree(r) == max_cover_degree}
            max_degree_cover = {p for r in max_degree_nodes for p in G[r]}
            G.remove_nodes_from(max_degree_nodes)
            G.remove_nodes_from(max_degree_cover)
            G.remove_nodes_from([r for r in G.prots() if G.degree(r) == 0])# remove unsupported
            cover.update(max_degree_nodes)
            for C in G.components():
                cover.update(greedy_minimal_cover(C))# recursive move
        return cover
    else:
        print(G)
        print(next(G.prots()))
        return {next(G.prots())} # the only protein left
