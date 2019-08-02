from .bigraph import BiGraph
from .min_cover import greedy_minimal_cover


class ProtPepGraph(BiGraph):
    def prots(self, deg=None):
        """Iterate proteins."""
        yield from self.A(deg)
    proteins = prots
    
    def prot_cnt(self, deg=None):
        return self.A_cnt(deg)

    def peps(self, deg=None):
        """Iterate peptides."""
        yield from self.B(deg)
    peptides = peps

    def pep_cnt(self, deg=None):
        return self.B_cnt(deg)

    def prot_pep_pairs(self):
        """Iterate potential peptide-protein explanations."""
        yield from self.AB()
    protein_peptide_pairs = prot_pep_pairs

    def __repr__(self):
        return "ProtPepGraph(proteins {} petpides {} links {})".format(*self.nodes_cnt(), len(self.edges))

    def remove_lonely_and_unsupported(self, min_pepNo_per_prot=2):
        """Get the full protein-peptide graph.
        
        Args:
            min_pepNo_per_prot (int): The minimal number of peptides per protein.
        Returns:
            Tuple: a tuple with lonely peptide and prots, and a tuple of peptides not belonging to unsupproted proteins and these proteins.    
        """
        # removing pairs r-p, where both r and p have no other neighbors
        lonely_peps = set(self.peps(deg=0))
        lonely_prots = set(self.prots(deg=0))
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

    def get_minimal_graph_2(self, max_iter=float('inf')):
        """Get the induced minimal graph.

        The minimal graph is a version of the PepProtGraph that merges nodes that share the same neighbors and reduces that further by finding the minimal number of protein groups needed to explain the observed peptide groups.

        Args:
            max_iter (int): The maximal number of repeats of the greedy minimal cover finder.

        Returns:
            tuple containing the induced minimal PepProtGraph and the set of protein groups not necessary to explain peptide groups.
        """
        G = self.form_groups()
        beckham_prot_groups = set([])
        prev_prots_cnt = G.prot_cnt()
        it = 0
        while it < max_iter:
            cover = greedy_minimal_cover(G)
            noncovering = {rg for rg in G.prots() if rg not in cover}
            beckham_prot_groups.update(noncovering)
            # removing beckham prots will not leave any peptides without neighbors, as they are covered by min_prot_cover.
            G.remove_nodes_from(noncovering)
            # after removing nodes from cover, some covered nodes might have lost proteins that discerned them from others.
            # Otherwise said, they will now have the same neighboring protein groups as some other peptide groups, which
            # was impossible before. To fix it, we have to reform the groups.
            G = G.form_groups(merging_merged=True)
            prots_cnt = G.prot_cnt()
            if prev_prots_cnt == prots_cnt:
                break
            else:
                prev_prots_cnt = prots_cnt
            it += 1
        return G, beckham_prot_groups

    def pop_unsupported(self):
        """Pop some of the unsupported proteins.

        Remove from the graph (and yield) proteins that do not have unique peptides
        and for any peptide they have there is a protein that has unique peptides.
        Later in the code, these are reffered to as inuprots surrounded by uniprots."""
        while True:
            supported = {r for p in self.peps(deg=1) for r in self[p]}
            shared_peps = {p for r in supported for p in self[r] if self.degree(p) > 1}
            if shared_peps:
                unsupported = {r for r in self.prots() if all(p in shared_peps for p in self[r])}
                if unsupported:
                    yield from unsupported
                    self.remove_nodes_from(unsupported)
                else:
                    break
            else:
                break

    def get_minimal_graph(self):
        """Get the induced minimal graph.

        The minimal graph is a version of the PepProtGraph that merges nodes that share the same neighbors and reduces that further by finding the minimal number of protein groups needed to explain the observed peptide groups.

        Returns:
            tuple containing the induced minimal PepProtGraph and the set of protein groups not necessary to explain peptide groups.
        """
        H = self.form_groups()
        cover = greedy_minimal_cover(H)
        noncovering = {r for r in H.prots() if r not in cover}
        H.remove_nodes_from(noncovering)
        H = H.form_groups(merging_merged=True)
        noncovering.update(H.pop_unsupported())
        H = H.form_groups(merging_merged=True)
        return H, noncovering


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
