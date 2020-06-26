import numpy as np
import pandas as pd
import networkx as nx
from networkx.algorithms import community
from networkx.algorithms.community.quality import modularity
from collections import defaultdict
from collections import Counter
from analysis_scripts import compute_analytics, graph_loader, projected_graph

def to_gephi(graph, out_path):
    with open(out_path, 'w') as f_out:
        f_out.write('graph water_molecules {\n')
        for node in list(graph.nodes(data=True)):
            f_out.write('%s [label=%s]\n' % (node[0], node[1]['label']))
        for e in list(graph.edges(data=True)):
            weight = float(e[2]['weight']) if 'weight' in e[2] else 1.0
            f_out.write('%s -- %s [label=%s, weight=%f]\n' % \
                    (e[0], e[1], e[2]['label'], weight))
        f_out.write('}\n')
    return
