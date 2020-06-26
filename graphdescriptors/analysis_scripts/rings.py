import numpy as np
import networkx as nx
import datetime
def get_ring_key(list_of_nodes):
    """
    Learn unique representations of a ring chain.
    :param list_of_nodes: a sequence of nodes that captures a ring
                          Example: [4, 1, 5, 4]
    :return a string representation of the path -> "1,4:1,5:4,5"
    """

    # Given a list of nodes V = [list_of_nodes[i-1], list_of_nodes[i]]
    # ','.join(sorted(list_of_nodes[i-1], list_of_nodes[i])) returns a
    # string of form "v1, v2" such that v1 = V_sorted[0] and v2  = V_sorted[1]

    list_of_ordered_edges = [','.join(sorted([str(list_of_nodes[i-1]), \
            str(list_of_nodes[i])])) \
            for i in range(1, len(list_of_nodes))]

    # For above input list list_of_ordered_edges = [(1,4), (1,5), (4,5)]
    # Next, sort the ordered edges and build a string of form "e1:e2:..:e_n" where
    # the e_i is encoded by the above statement.
    # ring_key for above example: "1,4:1,5:4,5"

    ring_key = ':'.join(sorted(list_of_ordered_edges))
    return ring_key

def find_rings(graph, ring_size, node=-1):
    """
    Find the number of rings associated with a node via depth-first-search.
    :param graph: Input graph to count rings in.
    :param node: Count number of rings attached to this node. Count rings in 
                 whole graph is node is set to -1.
    :param ring_size: Size of the ring, would be 3 for a trimer.
    :return Number of rings found
    """

    # Initialize visit stack.  We need a first-in-last-out data structure 
    # for depth-first-search (DFS). The stack will contain a list of paths that
    # will be expanded in a DFS fashion.
    visit_stack = [[n] for n in graph.nodes()] if node == -1 else [[node]]
    # Track visited nodes so we don't fall in a loop
    visited = set()
    rings = [] 
    ring_key_set = set()
   
    while len(visit_stack) > 0:
        path_to_expand = visit_stack.pop()
        #print('Path to expand: ')
        #print(path_to_expand)
        next_node = path_to_expand[-1]
        curr_path_len = len(path_to_expand)
        #print('Searching around node: %d' % next_node)

        for nbr in graph.neighbors(next_node):
            #print('Checking neighbor: %d' % nbr)
            if curr_path_len == ring_size:
                if nbr == path_to_expand[0]:
                    #print('We found a ring!')
                    ring_path = path_to_expand + [nbr]
                    ring_key = get_ring_key(ring_path)
                    G_sub=graph.subgraph(ring_path)
                    degrees=[row[1] for row in G_sub.degree()]

                    if 2 == np.mean(degrees):
                    # Check if this ring is already found.  
                    # For example, assumg we start the search from 
                    # node 0 in the test graph.  We should not double
                    # count the rings (0->1->3->2) and (0->2->3->1)
                        if ring_key not in ring_key_set:
                            rings.append(path_to_expand)
                            ring_key_set.add(ring_key)

            elif curr_path_len < ring_size:
                # Expand this path further provided there is no loop
                if nbr not in path_to_expand:
                    #print('     generating new path to expand:')
                    new_path_to_expand = path_to_expand + [nbr]
                    #print(new_path_to_expand)
                    visit_stack.append(new_path_to_expand)        

    # Upon termination 
    #print('Printing %d rings of size %d' % (len(rings), ring_size)) 
    #for r in rings:
    #    print(r)
    return rings

def enumerate_rings(graph):
    """
    Get the count of trimers, tetramers, pentamers and hexamers 
    """
    trimers = find_rings(graph, 3)
    tetramers = find_rings(graph, 4)
    pentamers = find_rings(graph, 5)
    hexamers = find_rings(graph, 6)
    return len(trimers), len(tetramers), len(pentamers), len(hexamers)
