# Greedy-Set-Cover(X, S) {
# U = X // U stores the uncovered items
# C = empty // C stores the sets of the cover
# while (U is nonempty) {
# select s[i] in S that covers the most elements of U
# add i to C
# remove the elements of s[i] from U
# }
# return C
# }
from pep2prot.graphs import random_bigraph

G = BiGraph([('a',0),('a',1),('b',1),('b',2),('c',2),('c',3),
	         ('d',3),('d',4),('e',4),('e',5),('f',6),('f',7),
	         ('g',7),('g',8),('g',9),('h',9),('h',10),('i',10),
	         ('i',11),('j',12),('j',14),('k',12),('k',13),
	         ('l',13),('l',14)])
G.draw()
# Facilitate the task by including the supported edges.
supported = BiGraph((a,b) for c in G.B() if G.degree(c)==1
                    for a in G[c] for b in G[a])
supported.draw(with_labels=True)
unsupported = BiGraph((a,b) for a,b in G.AB() if b not in supported)
# unsupported.draw(with_labels=True)
# input for MSC, if not a cycle

CC = next(unsupported.components())
CC.has_cycle()


TSC = BiGraph([('A',1),('A',2),('A',3),
              ('B',5),('B',2),('B',3),('B',4),
              ('C',4),('C',5),
              ('D',6),('D',7),('D',8),('D',9),
              ('E',8),('E',9),
              ('F',10),('F',12),
              ('G',11),('G',12),
              ('H',10),('H',11),('H',12)])
# G = random_bigraph(100, 50)
I.components()





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
	to_cover = G.B if A_covers_B else G.A
	covering = G.A if A_covers_B else G.B
	not_coverable = [m for m in to_cover() if G.degree(m)==0]
	G.remove_nodes_from(not_coverable)
	to_cover_cnt = G.B_cnt() if A_covers_B else G.A_cnt()
	cover = set([])
	while to_cover_cnt > 0:
		n_max = max(covering(), key=lambda n: len(G[n]))
		cover.add(n_max)
		n_max_neighbors = list(G[n_max])
		G.remove_nodes_from(n_max_neighbors)
		G.remove_node(n_max)
		G.remove_nodes_from([n for n in covering() if G.degree(n) == 0])#n_max covers what n does
		to_cover_cnt -= len(n_max_neighbors)
	return cover

TSC_MC = greedy_minimal_cover(TSC)
node_sizes = [40 if n in TSC_MC else 10 for n in TSC]
TSC.draw(node_size=node_sizes)
# min_cover = [greedy_minimal_cover(cc) for cc in R.components()]




# Add this a procedure to generate a random bipartite graph.
# pos = nx.bipartite_layout(G, G.A())
# nx.draw(G, pos)
# plt.show()


# G.draw()



%%time
try:
	X = nx.algorithms.find_cycle(G)
except nx.NetworkXNoCycle:
	pass


