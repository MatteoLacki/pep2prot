"""Class representing a bipartite graph."""

import networkx as nx
from collections import defaultdict
from networkx.algorithms import bipartite

from .min_set_cover import greedy_minimal_cover_2

set_union = lambda S: frozenset(r for s in S for r in s)

#TODO: it would be good to allow names in A and B to be the same???
class BiGraph(nx.Graph):
    """A bipartite graph."""
    def __init__(self, edges=None, *args, **kwds):
        """Instantiate a peptide-protein graph.

        Attention! Nodes in A should hash to something else than nodes in B.

        Arguments:
            edges (iterable): edges (a,b) s.t. a in A and b in B.
            *args, **kwds: arguments to the underlying networkx.Graph.
        """
        super().__init__(*args, **kwds)
        if edges:
            for a, z in edges:
                self.add_node(a, A=True)
                self.add_node(z, A=False)
                self.add_edge(a,z)

    def A(self):
        for n in self:
            if self.node[n]['A']:
                yield n

    def B(self):
        for n in self:
            if not self.node[n]['A']:
                yield n

    def add_A_node(self, node_for_adding, **attr):
        super().add_node(node_for_adding, A=True, **attr)

    def add_B_node(self, node_for_adding, **attr):
        super().add_node(node_for_adding, A=False, **attr)

    def add_A_nodes_from(self, nodes_for_adding, **attr):
        super().add_nodes_from(self, nodes_for_adding, A=True, **attr)

    def add_B_nodes_from(self, nodes_for_adding, **attr):
        super().add_nodes_from(self, nodes_for_adding, A=False, **attr)

    def add_AB_edge(self, a, b, **attr):
        if a in self:
            if not self.node[a]['A']:
                raise Exception("Node ({}) already in the B set!".format(a.__repr__()))
        else:
            self.add_node(a, A=True)
        if b in self:
            if self.node[b]['A']:
                raise Exception("Node ({}) already in the A set!".format(b.__repr__()))
        else:
            self.add_node(b, A=False)
        self.add_edge(a, b, **attr)

    def nodes_cnt(self):
        i = 0
        for a in self.A():
            i += 1
        return i, len(self) - i

    def A_cnt(self):
        i = 0
        for a in self.A():
            i += 1
        return i

    def B_cnt(self):
        i = 0
        for b in self.B():
            i += 1
        return i

    def __repr__(self):
        return "BiGraph(#A={} #B={} #E={})".format(*self.nodes_cnt(), len(self.edges))

    #TODO: start using bipartite layout for components.
    def draw(self, show=True, with_labels=False, node_size=10, *args, **kwds):
        import matplotlib.pyplot as plt
        node_colors = [('red' if self.node[n]['A'] else 'blue') for n in self]
        nx.draw_networkx(self,
                         node_color=node_colors,
                         with_labels=with_labels,
                         node_size=node_size,
                         *args, **kwds)
        if show:
            plt.show()

    def zero_degs(self):
        a0 = b0 = 0
        for n in self:
            if self.degree(n) == 0:
                if self.node[n]['A']:
                    a0 += 1
                else:
                    b0 += 1
        return a0, b0

    def cc_cnt(self):
        i = 0
        for _ in nx.connected_components(self):
            i += 1
        return i

    def AB(self):
        """Edge iteration that preserves the order of the reported node tuples."""
        for n,m in self.edges:
            if self.node[n]['A']:
                yield n,m
            else:
                yield m,n

    def components(self):
        for cc in nx.connected_components(self):
            yield self.subgraph(cc)

    def merge_nodes(self, AorB, merging_merged=False):
        """Merge nodes with the same neigbors in the other set.

        Calling with 'AorB == "A"', will merge 'A' nodes with common neighbors in 'B'.
        Calling with 'AorB == "B"', will merge 'B' nodes with common neighbors in 'A'.
        Merged nodes are kept as frozensets of nodes.
        If nodes were frozensets, it will try to merge them.
        So don't put frozensets as nodes, or all will explode.

        Args:
            AorB ('A' or 'B'): which set of nodes to merge?
            merging_merged (boolean): in case of second round of merging, don't create sets of sets, but just bigger sets of merged nodes.
        """
        assert AorB in ('A', 'B'), "A or B. Literaly."
        neighbors_of_nodes_to_merge = defaultdict(set)
        for n in (self.A() if AorB == 'A' else self.B()):
            neighbors_of_nodes_to_merge[frozenset(self[n])].add(n)
        agg = set_union if merging_merged else frozenset
        res = ((agg(merged),m) if AorB=='A' else (m,agg(merged))
               for to_merge, merged in neighbors_of_nodes_to_merge.items()
               for m in to_merge)
        return self.__class__(res)

    def form_groups(self, merging_merged=False):
        """Merge all nodes with common neighbors."""
        F = self.merge_nodes('B', merging_merged)
        return F.merge_nodes('A', merging_merged)

    def has_cycle(self):
        try:
            _ = nx.algorithms.find_cycle(self)
            return True
        except nx.NetworkXNoCycle:
            return False

    def greedy_minimal_cover(self, A_covers_B=True):
        return greedy_minimal_cover_2(self, A_covers_B)



def random_bigraph(maxA=20, maxB=40, prob=.05):
    """Generate a random BiGraph.

    Args:
        maxA (int): max number of nodes in A.
        maxB (int): max number of nodes in B.
        prob (float): probability of edge formation.
    Returns:
        BiGraph: A random BiGraph.
    """
    assert prob >= 0 and prob <= 1
    G = nx.algorithms.bipartite.generators.random_graph(maxA, maxB, p=prob)
    G.remove_nodes_from([n for n in G if G.degree(n) == 0])
    G = max((G.subgraph(cc) for cc in nx.connected_components(G)), key=lambda cc: len(cc))
    return BiGraph((a,b) for a in nx.bipartite.sets(G)[0] for b in G[a])



#TODO: modify the draw function to add a legend for colors.
#TODO: add the simplify method.
class ProtPepGraph(BiGraph):
    def prots(self):
        yield from self.A()

    def peps(self):
        yield from self.B()

    def prot_pep_pairs(self):
        yield from self.AB()

    def proteins(self):
        yield from self.prots()

    def peptides(self):
        yield from self.prots()

    def __repr__(self):
        return "ProtPepGraph(proteins {} petpides {} links {})".format(*self.nodes_cnt(), len(self.edges))
