"""Class representing a bipartite graph."""

import networkx as nx
from networkx.algorithms import bipartite


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

    def draw(self, show=True):
        import matplotlib.pyplot as plt
        node_colors = [('red' if self.node[n]['A'] else 'blue') for n in self]
        nx.draw_networkx(self, node_color=node_colors)
        if show:
            plt.show()

    def A(self):
        for n in self:
            if self.node[n]['A']:
                yield n

    def B(self):
        for n in self:
            if not self.node[n]['A']:
                yield n

    def draw(self, show=True):
        import matplotlib.pyplot as plt
        node_colors = [('red' if self.node[n]['A'] else 'blue') for n in self]
        nx.draw_networkx(self, node_color=node_colors)
        if show:
            plt.show()

    def nodes_cnt(self):
        i = 0
        for a in self.A():
            i += 1
        return i, len(self) - i
