import numpy as np
import pandas as pd
import networkx as nx

def read_lines(file_name):
    '''parse the output from the SchNet evalutation
       1) a list of the coordinates of each cluster (returns 1st)
          ['O   1.922131  -1.3179458  -2.891314\n', 'H   1.9396669  -1.812519  -2.0610204\n', ...]
       2) the actual energy (returns 2nd)
       3) the predicted energy (returns 3rd)'''
    
    # Load the output file and separate into energies, coordinates, and number of atoms
    coordinates_list=[]; energies_list=[]; atoms_list=[]
    file = open(file_name)
    
    for line in file:
        if 'Predicted' in line:
            energies_list.append(line)
        elif 'O' in line or 'H' in line:
            coordinates_list.append(line)
        else:
            atoms_list.append(int(line))        
    # Divide the structure list into clusters
    start=0; cluster_list=[]
    for i, atoms in enumerate(atoms_list):
        cluster_list.append([x for x in coordinates_list[start:start+atoms]])
        start += atoms
    
    return cluster_list, [float(str.split(x)[1]) for x in energies_list], [float(str.split(x)[3]) for x in energies_list]

def read_lines_base(file_name):
    '''parse the database xyz files
       this code only works if all structures have the same number of atoms
       1) a list of the coordinates of each cluster (returns 1st)
          ['O   1.922131  -1.3179458  -2.891314\n', 'H   1.9396669  -1.812519  -2.0610204\n', ...]
       2) the energy (returns 2nd)'''
    with open(file_name) as f:
        n_atoms = f.readline()       #atoms per cluster
        n_lines = 2 + int(n_atoms)   #lines per cluster (#atoms + energy + coords)
    with open(file_name) as f:
        lines = f.readlines()
    energies = np.array(lines[1::n_lines],dtype='float32') 
    structure_list, energy_list = [], []
    for n in range(int(energies.shape[0])):
        structure_list.append(lines[n_lines*(n+1)-n_lines:n_lines*(n+1)])
        structure_list[n][1]=float(structure_list[n][1])  #energy in float
        energy_list.append(float(structure_list[n][1]))
    return structure_list, energy_list


def adjacency_list(coord_list):
    '''Creates adjacency list to form graph in graph_analytics.
       Structure list must be ordered such that the two covalently
       bound hydrogens directly follow their oxygen.
       Definition of a hydrogen bond obtained from
       https://aip.scitation.org/doi/10.1063/1.2742385'''
    s = pd.Series(coord_list)
    dfx=s.str.split(expand=True)
    coords=[str.split(x) for x in coord_list]
    
    #delete atom label
    for j in coords:
        del j[0]

    cov_bonds, h_bonds, labels = [],[],[]
    for i, row in dfx.iterrows():
        labels.append([str(i+1),row[0]])
    q_1_2=[]
    for i in range(len(dfx)):
        if s[i].split()[0]=='O':
            cov_bonds.append([str(i+1),str(i+2),'covalent'])
            cov_bonds.append([str(i+1),str(i+3),'covalent'])
            h1=np.array(s[i+1].split()[1:],dtype='float64')
            h2=np.array(s[i+2].split()[1:],dtype='float64')
            o=np.array(s[i].split()[1:],dtype='float64')
            q_1_2.append([h1-o, h2-o])
    v_list=[np.cross(q1,q2) for (q1,q2) in q_1_2]
    for idx, v in enumerate(v_list):
        for index, both_roh in enumerate(q_1_2):
            for h_index, roh in enumerate(both_roh):
                indexO=((idx+1)*3)-2
                indexH=((index+1)*3)-2+(h_index+1)
                o_hbond = s[indexO-1].split()
                try:
                    h_hbond= s[indexH-1].split()  #will enumerate past the list if you let it
                except KeyError:
                    continue
                dist = np.linalg.norm(np.array(o_hbond[1:],dtype='float64')-np.array(h_hbond[1:],dtype='float64'))
                if (dist>1) & (dist<2.8):
                    angle = np.arccos(np.dot(roh, v)/(np.linalg.norm(roh)*np.linalg.norm(v)))*(180.0/np.pi)
                    if angle > 90.0:
                        angle=180.0-angle
                    N = np.exp(-np.linalg.norm(dist)/0.343)*(7.1-(0.05*angle)+(0.00021*(angle**2)))
                    if N >=0.0085:
                        h_bonds.append([str(indexO),str(indexH),'hydrogen'])

    return labels, cov_bonds, h_bonds, coords


def load_graph(struct):
    '''loads the graph for a single structure and returns the graph'''
    l,c,h,coords=adjacency_list(np.array(struct))
    node_labels = dict()
    for i in range(len(l)):
        node_labels[l[i][0]] = l[i][1]
    edges=c+h
    graph = nx.Graph()
    for k,v in node_labels.items():
        graph.add_node(k, label=v, coords=coords[int(k)-1])
    for triple in edges:
        atom1 = [float(x) for x in coords[int(triple[0])-1]]
        atom2 = [float(x) for x in coords[int(triple[1])-1]]        
        distance = np.linalg.norm(np.array(atom2)-np.array(atom1))
        graph.add_edge(triple[0], triple[1], label=triple[2], length=distance)
    return graph, node_labels, edges


