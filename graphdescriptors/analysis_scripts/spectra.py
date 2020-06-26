import numpy as np
import pandas as pd
import networkx as nx
from networkx.algorithms import community
from networkx.algorithms.community.quality import modularity
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
    return molecular_role_distribution
