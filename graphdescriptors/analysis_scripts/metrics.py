import numpy as np
import pandas as pd
import networkx as nx
from networkx.algorithms import community
from networkx.algorithms.community.quality import modularity
from collections import defaultdict
from collections import Counter
from analysis_scripts import compute_analytics, graph_loader, projected_graph


def compute_atoms(G):
    return len(G)

def compute_oxygens(G):
    projected_graph = projected_graph.project_oxygen_role_based_graph(G)
    return len(projected_graph)

def compute_hydrogens(G):
    projected_graph = projected_graph.project_oxygen_role_based_graph(G)
    return int(len(G)-len(projected_graph))

def compute_dangling_hydrogens(G):
    vals = [x for (y,x) in G.degree]
    degree1 = [x for x in vals if x==1]
    return degree1

def compute_aspl(G):
    '''If graph is not connected returns -1'''
    projected_graph = projected_graph.project_oxygen_role_based_graph(G)
    try:
       aspl= nx.average_shortest_path_length(projected_graph)
       return aspl
    except:
       return -1

def compute_wiener_index(G):
    '''If graph is not connected returns -1'''
    projected_graph = projected_graph.project_oxygen_role_based_graph(G)
    try:
       wiener = nx.wiener_index(projected_graph)
       return wiener
    except:
       return -1

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
