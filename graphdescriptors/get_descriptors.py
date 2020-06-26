import numpy as np
import pandas as pd
import networkx as nx
from networkx.algorithms import community
from networkx.algorithms.community.quality import modularity
from collections import defaultdict
from collections import Counter
from analysis_scripts import compute_analytics, graph_loader, projected_graph, rings, metrics
import sys
import argparse
import os.path as op

p = argparse.ArgumentParser()
p.add_argument("--data_path", type=str, required=True, help="Path to xyz file to analyze.")
p.add_argument("--output", type=str, required=True, help="Path + file name for output df (end in csv).")
p.add_argument("--min_dir", type=str, required=True, help="Path to xyz files of lowest energy structures.")
args = p.parse_args()

cluster_list, actual_energy, predicted_energy = graph_loader.read_lines(args.data_path)
print(str(len(cluster_list))+' clusters found')

a = [compute_analytics.graph_analytics(x) for x in cluster_list]
a = np.swapaxes(a, 0, 1)

similarity_list, projected_similarity_list=[],[]
trimers,tetramers,pentamers,hexamers=[],[],[],[]
degrees=[]
for i, cluster in enumerate(cluster_list):
    cluster_size = int(int(len(cluster))/3)
    xyzfile= f'W{cluster_size}_geoms_all_lowestE.xyz' 
    base_cluster_list,_ = graph_loader.read_lines_base(op.join(args.min_dir,xyzfile))
    #load lowest energy structure from the full list
    G_base,_,_ = graph_loader.load_graph(base_cluster_list[0][2:]) 
    proj_G_base = projected_graph.project_oxygen_role_based_graph(G_base)
    try:
       #Similarity of full graph
       G_variable,_,_ = graph_loader.load_graph(cluster)
       similarity_list.append(compute_analytics.calculate_similarity(G_base, G_variable))
       #Similarity of oxygen projected graph
       proj_G_variable = projected_graph.project_oxygen_role_based_graph(G_variable)
       projected_similarity_list.append(compute_analytics.calculate_similarity(proj_G_base, proj_G_variable))
       tri, tet, pent, hexa = rings.enumerate_rings(proj_G_variable)
       trimers.append(tri)
       tetramers.append(tet)
       pentamers.append(pent)
       hexamers.append(hexa)
       degs=[v for k,v in proj_G_variable.degree()]
       degrees.append(np.mean(degs))
    except:
       #disconnected graphs fail test and get s values of -1
       projected_similarity_list.append(-1)
       similarity_list.append(-1)


d = {'Nodes': a[0], 'Edges': a[1], 'Diameter': a[2],
     'Dangling Hydrogens': a[3], 'Average Shortest Path Length': a[4],
     'daa': a[5], 'dda': a[6], 'da': a[7], 'aa': a[8], 'dd': a[9], 
     'ddaaa': a[10], 'ddaa': a[11], 'a': a[12], 'd': a[13],
     'Wiener Index': a[14], 'N Communities': a[15], 'Modularity': a[16], 
     'Actual Energy': actual_energy,'Predicted Energy': predicted_energy,
     'Similarity': similarity_list,'Projected Similarity': projected_similarity_list,
     'Trimers': trimers, 'Tetramers': tetramers, 'Pentamers': pentamers,
     'Hexamers': hexamers, 'Degree': degrees}

df = pd.DataFrame(d)

df['Cluster Size'] = df['Nodes']/3
df['H-bonds'] = df['Edges']-(df['Cluster Size']*2)
df['H-bond Percent'] = df['H-bonds']/df['Edges']
df['Unsaturation'] = df['H-bonds']/df['Dangling Hydrogens']
df['daa Percent'] = df['daa']/df['Nodes']
df['dda Percent'] = df['dda']/df['Nodes']
df['Absolute Error per Water'] = ((df['Actual Energy']-df['Predicted Energy'])/df['Cluster Size']).apply(lambda x: np.abs(x))
df['Percent Error']=((df['Actual Energy']-df['Predicted Energy'])/df['Actual Energy']).apply(lambda x: 100*np.abs(x))
df['H-bond Saturation']=(df['H-bonds']/df['Cluster Size'].apply(lambda x: x*4)).apply(lambda x: x*2)

df.to_csv(args.output, index=False)
