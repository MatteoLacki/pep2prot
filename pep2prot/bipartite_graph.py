"""Class representing a bipartite graph."""

import networkx as nx
from collections import defaultdict
# from networkx.algorithms import bipartite


class BiGraph(nx.Graph):
    """A bipartite graph."""
    def __init__(self, edges=None, *args, **kwds):
        """Instantiate a peptide-protein graph.

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

    def __repr__(self):
        return "BiGraph(#A={} #B={} #E={})".format(*self.nodes_cnt(), len(self.edges))

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

    def merge_nodes(self, AorB):
        """Merge nodes with the same neigbors in the other set.

        Calling with 'AorB == "A"', will merge 'A' nodes with common neighbors in 'B'.
        Calling with 'AorB == "B"', will merge 'B' nodes with common neighbors in 'A'.
        """
        assert AorB in ('A', 'B'), "A or B. Literaly."
        neighbors_of_nodes_to_merge = defaultdict(set)
        for n in (self.A() if AorB == 'A' else self.B()):
            neighbors_of_nodes_to_merge[frozenset(self[n])].add(n)
        return self.__class__(
            (frozenset(merged),m) if AorB == 'A' else (m,frozenset(merged))
            for to_merge, merged in neighbors_of_nodes_to_merge.items()
            for m in to_merge
        )

    def form_groups(self):
        """Merge all nodes with common neighbors."""
        F = self.merge_nodes('B')
        return F.merge_nodes('A')


class ProtPepGraph(BiGraph):
    def prots(self):
        yield from self.A()

    def peps(self):
        yield from self.B()

    def prot_pep_pairs(self):
        yield from self.AB()

    def __repr__(self):
        return "ProtPepGraph(#A={} #B={} #E={})".format(*self.nodes_cnt(), len(self.edges))
