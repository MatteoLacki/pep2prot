"""Class representing a bipartite graph."""

import networkx as nx
from collections import defaultdict
from networkx.algorithms import bipartite

from .min_set_cover import greedy_minimal_cover_2

set_union = lambda S: frozenset(r for s in S for r in s)


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
            merging_merged (boolean): in case of second round of merging, don't create sets of sets, but just bigger sets of merged nodes. If the merged nodes were {a,b}, {c,d}, then don't create {{a,b}, {c,d}}, but rather {a,b,c,d}.
        """
        assert AorB in ('A', 'B'), "A or B. Literaly."
        # merge nodes that share the same neighbors
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

    def minimal_cover(self, A_covers_B=True):
        return greedy_minimal_cover_2(self, A_covers_B)

    @classmethod
    def random(cls, maxA=20, maxB=40, prob=.05):
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
        return cls((a,b) for a in nx.bipartite.sets(G)[0] for b in G[a])



class ProtPepGraph(BiGraph):
    def prots(self):
        """Iterate proteins."""
        yield from self.A()

    def peps(self):
        """Iterate peptides."""
        yield from self.B()

    def prot_pep_pairs(self):
        """Iterate potential peptide-protein explanations."""
        yield from self.AB()

    def proteins(self):
        """Iterate proteins."""
        yield from self.prots()

    def peptides(self):
        """Iterate peptides."""
        yield from self.prots()

    def __repr__(self):
        return "ProtPepGraph(proteins {} petpides {} links {})".format(*self.nodes_cnt(), len(self.edges))

    def minimal_protein_cover(self):
        return self.minimal_cover(A_covers_B=True) # A = proteins, B = peptides, as seen above...

    def remove_lonely_and_unsupported(self, min_pepNo_per_prot=2):
        """Get the full protein-peptide graph.
        
        Args:
            min_pepNo_per_prot (int): The minimal number of peptides per protein.
        Returns:
            Tuple: a tuple with lonely peptide and prots, and a tuple of peptides not belonging to unsupproted proteins and these proteins.    
        """
        # removing pairs r-p, where both r and p have no other neighbors
        lonely_peps = {p for p in self.peps() if self.degree(p) == 0}
        lonely_prots = {p for p in self.prots() if self.degree(p) == 0}
        self.remove_nodes_from(lonely_peps)
        self.remove_nodes_from(lonely_prots)

        prots_no_peps = {r for r in self.prots() if self.degree(r) < min_pepNo_per_prot}
        peps_no_prots = {p for r in prots_no_peps for p in self[r]}
        self.remove_nodes_from(prots_no_peps)
        # Remove only those peps that did not have any other neighbors than low count prots.
        # e.g. A-0-B C-1-D-2 and minimal number of peptide per protein is 2, then prots_no_peps = {A,B,C} and peps_no_prots = {0}
        peps_no_prots = {p for p in peps_no_prots if self.degree(p) == 0}
        self.remove_nodes_from(peps_no_prots)
        return (lonely_peps, lonely_prots), (peps_no_prots, prots_no_peps)

    def get_minimal_graph(self):
        """Get the induced minimal graph.

        The minimal graph is a version of the PepProtGraph that merges nodes that share the same neighbors and reduces that further by finding the minimal number of protein groups needed to explain the observed peptide groups.

        Returns:
            tuple containing the induced minimal PepProtGraph and the set of protein groups not necessary to explain peptide groups.
        """
        min_self = self.form_groups()
        min_prot_cover = min_self.minimal_protein_cover()
        beckham_prot_groups = {rg for rg in min_self.prots() if rg not in min_prot_cover}
        # removing beckham prots will not leave any peptides without neighbors, as they are covered by min_prot_cover.
        min_self.remove_nodes_from(beckham_prot_groups)
        # after removing protein groups some peptide groups might have lost proteins that discerned them from others.
        # Otherwise said, they will now have the same neighboring protein groups as some other peptide groups, which
        # was impossible before. To fix it, we have to reform the groups.
        min_self = min_self.form_groups(merging_merged=True)
        return min_self, beckham_prot_groups


# TODO: update this test.
# def test_node_merger():
#     """TESTING IF ALL NODES HAVE DIFFERENT NEIGHBORS."""
#     min_pepNo_per_prot = 2
#     max_rt_deviation = 1
#     path = Path(r"~/Projects/pep2prot/pep2prot/data").expanduser()
#     D = read_isoquant_peptide_report(path/'peptide_report.csv')
#     D, I_cols = preprocess_isoquant_peptide_report(D)
#     unique_columns = ['peptide_overall_max_score','peptide_fdr_level','peptide_overall_replication_rate','prots','pre_homology_accessions','pi','mw']
#     DD = complex_cluster_buster(D, I_cols, unique_columns, max_rt_deviation)
#     prot2seq = {r for rg in D2.prots for r in rg}
#     prot2seq = read_n_check_fastas(path/'mouse.fasta', prot2seq)
#     H, RWEP, BRR = get_peptide_protein_graph(DD)
#     x = Counter(frozenset(H[r]) for r in H.prots())
#     assert set(Counter(x.values())) == {1}
#     x = Counter(frozenset(H[r]) for r in H.peps())
#     assert set(Counter(x.values())) == {1}
