import numpy as np
import pandas as pd
import networkx as nx
from collections import defaultdict
from collections import Counter
from analysis_scripts import compute_analytics, graph_loader, projected_graph

def graph_projection(graph):
    """
    Project a new graph based on connectivity between a subset of nodes in the
    input graph (subsequently called "targets") via another subset (non-overlapping
    with the targets) of nodes referred to as "connectors". Connectors can be
    the complementary set of targets, aka set(V(graph) - targets) or a subset of that
    complementary set.
    """
    targets = [node for node,attr in graph.nodes(data=True) if attr['label']=='O']
    connectors = set([node for node,attr in graph.nodes(data=True) if attr['label']=='H'])
    weighted_edges = Counter()

    for node1 in connectors:
        neighbors = graph[node1]
        neighboring_targets = []
        for node2 in neighbors:
            if node2 in targets:
                neighboring_targets.append(node2)
        # iterate over the list of all candidate nodes that are connected to
        node_count = len(neighboring_targets)
        for i in range(node_count-1):
            for j in range(i+1, node_count):
                weighted_edges[(neighboring_targets[i], neighboring_targets[j])] += 1
    return weighted_edges



def project_oxygen_role_based_graph(graph):
    """
    Computes a graph capturing the interaction between oxygen atoms in a network
    of water molecules.  Nodes in the output graph represent oxygen atoms, and
    are labeled based on their participations in covalent and hydrogen bonds.
    An oxygen atom that is connected to other hydrogen atoms via two covalent
    bond and one hydrogen bond would be labeld as "c2:h1".
    """
    node_roles, role_distribution = compute_analytics.molecular_role_distribution(graph)
    weighted_edges = graph_projection(graph)

    projected_graph = nx.Graph()
    for k,v in node_roles.items():
        projected_graph.add_node(k, label=v)
    for edge, count in weighted_edges.items():
        projected_graph.add_edge(edge[0], edge[1], label='hydrogen bond')

    return projected_graph#compute_analytics.projected_oxygen_graph_metrics(projected_graph)
