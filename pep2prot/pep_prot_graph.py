"""Class representing the petpide-protein graph."""

import networkx as nx
from collections import Counter

# input for the data processing? Or maybe we should work on multiindices????
# no, because further on we will use columns and it's inconvenient. Use strings
# assumption: prots contains a column with proteins separated by some string.
# could very well modify the format for broader generality.


class PepProtGraph(nx.Graph):
    """Extension of the networkx.Graph for peptide-protein inference.

    Internally, peptides are strings beginning with a dollar.
    No, not Yuan. A dollar.
    """
    def __init__(self, peptide2proteins=None, *args, **kwds):
        """Instantiate a peptide-protein graph.

        Arguments:
            peptide2proteins (iterable): map between petides and list of corresponding proteins.
            *args, **kwds: arguments to the underlying networkx.Graph.
        """
        super().__init__(*args, **kwds)
        if peptide2proteins:
            self.add_edges_from(("$"+p,r) for p, prots in peptide2proteins for r in prots)

    def prots(self):
        for n in self:
            if n[0] != "$":
                yield n

    def peps(self):
        for n in self:
            if n[0] == "$":
                yield n[1:]
    
    def _peps(self):
        for n in self:
            if n[0] == "$":
                yield n

    def is_prot(self, n):
        return n[0] != "$"

    def is_pep(self, n):
        return n[0] == "$"

    def peps_prots_cnt(self):
        p = r = 0
        for n in self:
            if self.is_pep(n):
                p += 1
            else:
                r += 1
        return p, r

    def low_deg_nodes(self, max_deg, node_selector):
        for n in self:
            if node_selector(n) and self.degree(n) <= max_deg:
                yield n

    def __repr__(self):
        return "pep_{}_prot_{}_graph: N={} E={}".format(*self.peps_prots_cnt(), len(self.nodes), len(self.edges))

    def del_prots_with_little_peps(self, min_peps_per_prot):
        self.remove_nodes_from(list(self.low_deg_nodes(min_peps_per_prot-1, self.is_prot)))
        self.remove_nodes_from(list(self.low_deg_nodes(0, self.is_pep)))

    def draw(self, show=True):
        import matplotlib.pyplot as plt
        node_colors = [('red' if self.is_prot(n) else 'blue') for n in self]
        nx.draw_networkx(self, node_color=node_colors)
        if show:
            plt.show()

    def pep_degs(self):
        return Counter(self.degree(p) for p in self._peps())