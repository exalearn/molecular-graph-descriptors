import numpy as np
import pandas as pd
import networkx as nx
from networkx.algorithms import community
from networkx.algorithms.community.quality import modularity
from collections import defaultdict
from collections import Counter
from analysis_scripts import compute_analytics, graph_loader, projected_graph

def molecular_role_distribution(graph):
    """
    Function that computes a donor-acceptor profile for each water molecule.
    It accepts a graph with hydrogen and oxygen atoms as nodes, and returns
    two dictionaries.
    The first dictionary provides molecule-level role description. Keys in
    the dictionary correspond to the oxygen atom ids in each water molecule.
    The values are the corresponding donor-acceptor roles aggregated from the
    hydrogen and oxygen atoms in that molecule.  Profiles such as
    donor-acceptor-donor or donor-donor-acceptor are both represented as a1d2.
    The second dictionary is a Counter that aggregates the counts from
    the distribution stored in the first dictionary.
    Example: Counter({'a2d2': 12, 'a1d2': 9, 'a2d1': 8, 'a3d2': 1})
    """
    atomic_roles = defaultdict(Counter)
    node_labels = nx.get_node_attributes(graph, 'label')
    # Iterate over all hydrogen bonds.  Assign the hydrogen atom a "donor"
    # role and the oxygen atom an "acceptor" role.
    for edge in list(graph.edges(data=True)):
        if edge[2]['label'] == 'hydrogen':
            role = 'd' if node_labels[edge[0]] == 'H' else 'a'
            atomic_roles[edge[0]][role] += 1
            role = 'd' if node_labels[edge[1]] == 'H' else 'a'
            atomic_roles[edge[1]][role] += 1

    # Iterate over all oxygen atoms. For each oxygen atom, aggregate the
    # role counts from each of the hydrogen atoms connected to it via a
    # covalent bond. Thus each oxygen atom aggregates the roles for the
    # water molecule it is part of.
    molecular_roles = dict()
    molecular_role_distribution = Counter()
    for node in list(graph.nodes(data=True)):
        # Iterate over a list such as
        # [('24', {'label': 'H'}), ('25', {'label': 'O'})
        # Only do the computation of nodes representing oxygen atoms
        node_id = node[0]
        node_label = node[1]['label']
        if node_label != 'O':
            continue
        neighbors = graph[node_id]
        for n in neighbors:
            if neighbors[n]['label'] == 'covalent':
                roles = atomic_roles[n]
                atomic_roles[node_id] += roles
        tmp_list = []
        for k in sorted(atomic_roles[node_id].keys()):
            tmp_list.append('%s%d' % (k, atomic_roles[node_id][k]))
        role_key = ''.join(tmp_list)
        molecular_roles[node_id] = role_key
        molecular_role_distribution[role_key] += 1
    return molecular_roles,molecular_role_distribution


def projected_oxygen_graph_metrics(projected_graph):
    communities = community.greedy_modularity_communities(projected_graph)
    mod_score = modularity(projected_graph, communities) #not returned
    try:
       aspl = nx.average_shortest_path_length(projected_graph)
       wiener = nx.wiener_index(projected_graph)
    except:
       aspl = -1
       wiener = -1
    return aspl, wiener, len(communities), communities, mod_score



def graph_analytics(struct):
    '''loads the graph for a single structure and returns metrics in the following order
       1) number of nodes
       2) number of edges
       3) diameter
       4) number of nodes of degree 1
       5) average shortest path length -- O PROJECTED GRAPH
       6-12) number of donor/acceptor types -- Os only
       13) Wiener index -- O PROJECTED GRAPH
       14-15) number of communities and communities
    '''
    graph, node_labels, edges = graph_loader.load_graph(struct)
    vals = [x for (y,x) in graph.degree]
    degree1 = [x for x in vals if x==1]
    _, mrd_dict = molecular_role_distribution(graph)
    # Try/except in the case of an unconnected graph
    try:
        dia = nx.diameter(graph)
    except:
        dia = -1000
    daa, dda, da, aa, dd, ddaaa, ddaa, a, d = mrd_dict['a2d1'], mrd_dict['a1d2'], mrd_dict['a1d1'], mrd_dict['a2'], mrd_dict['d2'], mrd_dict['a3d2'], mrd_dict['a2d2'], mrd_dict['a1'], mrd_dict['d1']
    avg_short_path_length, w_index_O, n_communities, comms, mod_score = projected_oxygen_graph_metrics(projected_graph.project_oxygen_role_based_graph(graph))
    return len(node_labels), len(edges), dia, len(degree1), avg_short_path_length, daa, dda, da, aa, dd, ddaaa, ddaa, a, d, w_index_O, n_communities, mod_score#, comms

'''
   Eigenvector Similarity. Calculate the Laplacian eigenvalues for the adjacency matrices of each graph.
   For each graph, find the smallest k such that the sum of the k largest eigenvalues constitutes at least
   90% of the sum of all of the eigenvalues. If the values of k are different between the two graphs, then
   use the smaller one. The similarity metric is then the sum of the squared differences between the largest
   k eigenvalues between the graphs. This will produce a similarity metric in the range [0, ~H~^), where values
   closer to zero are more similar.
   https://stackoverflow.com/questions/12122021/python-implementation-of-a-graph-similarity-grading-algorithm
   https://www.cs.cmu.edu/~jingx/docs/DBreport.pdf
'''

def select_k(spectrum, minimum_energy = 0.9):
    running_total = 0.0
    total = sum(spectrum)
    if total == 0.0:
        return len(spectrum)
    for i in range(len(spectrum)):
        running_total += spectrum[i]
        if running_total / total >= minimum_energy:
            return i + 1
    return len(spectrum)

def calculate_similarity(G_base, G_variable):
    laplacian1 = nx.spectrum.laplacian_spectrum(G_base)
    laplacian2 = nx.spectrum.laplacian_spectrum(G_variable)

    k1 = select_k(laplacian1)
    k2 = select_k(laplacian2)
    k = min(k1, k2)

    similarity = sum((laplacian1[:k] - laplacian2[:k])**2)
    return similarity
