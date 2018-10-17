import sys
import os
import numpy as np
import networkx as nx
sys.path.append(os.path.realpath(__file__).rsplit('/', 2)[0])
import biographs as bg
import matplotlib.pyplot as plt
import amino_acids_conversion
import pandas as pd
from collections import deque
from matplotlib.colors import Normalize


def create_network(pdb_id, database=None, net=None, pdbs_path=None,
                   database_path=None, pathogenic=[], non_pathogenic=[],
                   both=[], colors='mutations',
                   sizes='neighborhood_watch_sharp'):
    
    if pdbs_path and database_path:
        if 'pdb' + pdb_id + '.pdb' in os.listdir(pdbs_path):
            pdb = os.path.join(pdbs_path, 'pdb' + pdb_id + '.pdb')
        else:
            pdb = os.path.join(pdbs_path, 'pdb' + pdb_id + '.ent')
        
        mol = bg.Pmolecule(pdb)
        net = mol.network()
        
        database = pd.DataFrame(pd.read_csv(database_path))
    
    node_labels = {}
    for node in net.nodes:
        info = database[database['Residue name'] == 'pdb' + pdb_id + node]
        if len(info) > 1: info = info.iloc[0] # check why more than one
        type_aa = amino_acids_conversion.three2one(info["Type of residue"].item())
        label = type_aa + node[1::] + ":" + node[0]
        node_labels[node] = label


    mutation_type = []
    neighborhood_watch_sharp = []
    neighborhood_watch_smooth = []
    degree = []
    pairwise_sharp = []
        
    for node in net.nodes:
        if colors == 'mutations' or sizes == 'mutations':
            seq_pos = node[1::]
            if seq_pos in pathogenic:
                mutation_type.append('tomato')
            elif seq_pos in non_pathogenic:
                mutation_type.append('limegreen')
            elif seq_pos in both:
                mutation_type.append('gold')
            else:
                mutation_type.append('lightgrey')
        
        elif colors or sizes:            
            k = nx.degree(net, node)            
            degree.append(k)
            weight = nx.degree(net, node, weight='weight')
            if colors == 'neighborhood_watch_sharp':
                if weight / k < 5: neighborhood_watch_sharp.append('blue')
                elif weight / k < 10:  neighborhood_watch_sharp.append('cyan')
                elif weight / k < 15:  neighborhood_watch_sharp.append('greenyellow')
                elif weight / k < 20:  neighborhood_watch_sharp.append('yellow')
                elif weight / k < 25:  neighborhood_watch_sharp.append('orange')
                else: neighborhood_watch_sharp.append('red')
            
            elif sizes == 'neighborhood_watch_sharp':
                if weight / k < 5: neighborhood_watch_sharp.append(1000)
                elif weight / k < 10:  neighborhood_watch_sharp.append(4000)
                elif weight / k < 15:  neighborhood_watch_sharp.append(7000)
                elif weight / k < 20:  neighborhood_watch_sharp.append(10000)
                elif weight / k < 25:  neighborhood_watch_sharp.append(13000)
                else: neighborhood_watch_sharp.append(16000)
            
            elif colors == 'neighborhood_watch_smooth' or sizes == 'neighborhood_watch_smooth':
                neighborhood_watch_smooth.append(weight / k)
    
    if colors == 'pairwise_sharp' or sizes == 'pairwise_sharp':
        for u, v in net.edges:
            wij = net.get_edge_data(u, v)['weight']
            if wij < 10: pairwise_sharp.append('blue')
            elif wij < 20:  pairwise_sharp.append('cyan')
            elif wij < 30:  pairwise_sharp.append('greenyellow')
            elif wij < 40:  pairwise_sharp.append('yellow')
            elif wij < 50:  pairwise_sharp.append('orange')
            else: pairwise_sharp.append('red')
                    
    color_map = []
    edge_color_map = []
    if colors == 'mutations': color_map = mutation_type
    elif colors == 'degree': color_map = degree
    elif colors == 'neighborhood_watch_sharp': color_map = neighborhood_watch_sharp
    elif colors == 'neighborhood_watch_smooth': color_map = neighborhood_watch_smooth 
    elif colors == 'pairwise_sharp': edge_color_map = pairwise_sharp
    
    size_map = []
#    edge_size_map = []
    if sizes == 'mutations': size_map = mutation_type
    elif sizes == 'degree': size_map = degree
    elif sizes == 'neighborhood_watch_sharp': size_map = neighborhood_watch_sharp
    elif sizes == 'neighborhood_watch_smooth': size_map = neighborhood_watch_smooth
#    elif sizes == 'pairwise_sharp': edge_size_map = pairwise_sharp
    
    return net, node_labels, size_map, color_map, edge_color_map




def draw_network(net, node_labels, sizes=None, color_map=None, cmap=None,
                 edge_color_map=None, name='', draw_edges=True,
                 draw_labels=True, pos=None, layout='spring',
                 default_size=2000, figsize=(130, 100), facecolor='w',
                 node_edge_colors=None, labels_size=30, node_edge_size=4,
                 pos_label=None):
    fig = plt.figure(figsize=figsize, facecolor=facecolor)
#    fig = plt.figure(facecolor=facecolor)
    
    if not pos:
        if layout == 'circular':
            pos=nx.circular_layout(net)
        else: # TO DO: put all cases 
            pos=nx.spring_layout(net)
    width = nx.get_edge_attributes(net, 'weight')
    
    if not sizes:
        sizes = default_size
    if not color_map:
        color_map = 'w'
    if not edge_color_map:
        edge_color_map = 'grey'
    
    if not node_edge_colors:
        node_edge_colors = ['k' for node in net.nodes]
        
    nx.draw_networkx_nodes(net, pos, node_size=sizes, node_color=color_map,
                           linewidths=node_edge_size, cmap=cmap,
                           edgecolors=node_edge_colors)
    if draw_labels:
        if not pos_label:
            pos_label = pos
        nx.draw_networkx_labels(net, pos_label, labels=node_labels,
                                font_size=labels_size, font_weight='bold')
    if draw_edges:
        nx.draw_networkx_edges(net, pos, width=[w / 2 for w in width.values()],
                                                edge_color=edge_color_map)
    
    plt.axis('off')
#    plt.tight_layout()
#    fig.patch.set_facecolor(facecolor)
    plt.axis('equal')
    fig.savefig(name + '.pdf', facecolor=facecolor, bbox_inches='tight')
    plt.close()
    
def type_structure(node, neighbor):
    node_pos = int(node[1::])
    node_chain = node[0]
    neighbor_pos = int(neighbor[1::])
    neighbor_chain = neighbor[0]
    dist = neighbor_pos - node_pos
    abs_dist = abs(dist)
    if node_chain != neighbor_chain:
        structure = '4D'
        dist = 99999
    elif abs_dist == 1:
        structure = '1D'
    elif abs_dist <= 4:
        structure = '2D'
    else:
        structure = '3D'
    return (structure, dist)

def create_ego_network(net, central_node, all_labels, n_cases, sort,
                       pathogenic=[],
                       non_pathogenic=[], both=[], threshold_edge=20):
        ego = nx.ego_graph(net, central_node)
        color_map = []
        sizes = []
        labels = {}
        for node in ego.nodes:
            seq_pos = node[1::]
            if seq_pos in pathogenic:
                color_map.append('tomato')
            elif seq_pos in non_pathogenic:
                color_map.append('limegreen')
            elif seq_pos in both:
                color_map.append('gold')
            else:
                color_map.append('lightgrey')
             
            k = nx.degree(net, node)
            weight = nx.degree(net, node, weight='weight')
            if weight / k < 10: sizes.append(200)
            elif weight / k > 16: sizes.append(1000)
            else: sizes.append(3500)
            labels[node] = all_labels[node]
        neighbors = ego.copy()
        neighbors.remove_node(central_node)
        neighbors_structure = {}
        node_borders = {}
        if n_cases and len(n_cases) > 0:
            check_colors = True
            edge_colors = {}
        else:
            check_colors = False
            edge_colors = 'grey'

        for neighbor in neighbors:
            struct, dist = type_structure(central_node, neighbor) # to modify
            neighbors_structure[neighbor] = (struct, dist)
            if struct == '1D': node_borders[neighbor] = 'turquoise'
            elif struct == '2D': node_borders[neighbor] = 'turquoise'
            elif struct == '3D': node_borders[neighbor] = 'cadetblue'
            else: node_borders[neighbor] = 'cadetblue'
            
            if check_colors:
                num = n_cases[neighbor]
                if num < threshold_edge:
                    edge_colors[(central_node, neighbor)] = 'red'
                else:
                    edge_colors[(central_node, neighbor)] = 'green'
        
        if sort == True:
            sorted_neighbors = deque(sorted(neighbors_structure,
                                            key=lambda x: neighbors_structure[x][1]))
            initial = central_node[0] + str(int(central_node[1::]) + 1)
            if initial in sorted_neighbors:
                sorted_neighbors.rotate(-sorted_neighbors.index(initial))
            else:
                initial = central_node[0] + str(int(central_node[1::]) - 1)
                sorted_neighbors.rotate(-sorted_neighbors.index(initial) + 1)
        
        else:
            sorted_neighbors = list(neighbors)
        
        pos_original=nx.circular_layout(neighbors)
        pos = {}
        for i in range(len(sorted_neighbors)):
            pos[sorted_neighbors[i]] = list(pos_original.values())[i]
        pos[central_node] = np.array([0, 0])
        node_borders[central_node] = 'grey'
        width = nx.get_edge_attributes(ego, 'weight')
        
        if check_colors:
            edges = ego.edges()
            edge_colors_list = []
            for u, v in edges:
                try:
                    edge_colors_list.append(edge_colors[(u, v)])
                except KeyError:
                    try:
                        edge_colors_list.append(edge_colors[(v, u)])
                    except KeyError:
                        edge_colors_list.append('grey')
            edge_colors = edge_colors_list
        
        return ego, labels, pos, sizes, width, color_map, node_borders, edge_colors

#
#def draw_many_neighborhoods(net, sequence, all_labels, start=0, stop=None, pathogenic=[], non_pathogenic=[], both=[], save_fig=True, file_name='figure', n_cases=None, threshold_edge=20):
#    index = 0
#    plt.figure(figsize=(15*1.5, 10*1.5))
#    if stop == None: stop = len(sequence) 
#    for pos in sequence[start:stop]:
#        try:
#            int(pos[0])
#            central_node = 'A' + pos
#        except ValueError:
#            central_node = pos
#        if central_node in net.nodes:
#            ego, labels, pos, sizes, width, color_map, node_borders, edge_colors = create_ego_network(net, central_node, all_labels, pathogenic, non_pathogenic, both, n_cases, threshold_edge)
#            index += 1
#            plt.subplot(3, 4, index)
#            nx.draw_networkx_nodes(ego, pos, node_size=sizes, node_color=color_map, edgecolors=[node_borders[n] for n in ego.nodes], linewidths=4)
#            nx.draw_networkx_labels(ego, pos, labels=labels, font_size=12, font_weight='bold')
#            nx.draw_networkx_edges(ego, pos, width=[w / 5 for w in width.values()], edge_color=edge_colors)
#    #        nx.draw_networkx_edge_labels(ego, pos, node_size=sizes, node_color=color_map, edge_labels=width)
#            plt.axis('off')
#    plt.tight_layout()
#    plt.savefig(file_name + ".pdf")
#    return ego, labels, pos, sizes, width, color_map, node_borders

def draw_neighborhood(net, pos, all_labels, pathogenic=[],
                      non_pathogenic=[], both=[], save_fig=True,
                      file_name='figure', n_cases=None, threshold_edge=20,
                      sort=True):
    plt.figure(figsize=(6, 6))
    try:
        int(pos[0])
        central_node = 'A' + pos
    except ValueError:
        central_node = pos
    if central_node in net.nodes:
        ego, labels, pos, sizes, width, color_map, node_borders, edge_colors = create_ego_network(net, central_node, all_labels, n_cases, pathogenic, non_pathogenic, both, threshold_edge=threshold_edge)
        nx.draw_networkx_nodes(ego, pos, node_size=sizes, node_color=color_map, edgecolors=[node_borders[n] for n in ego.nodes], linewidths=4)
        nx.draw_networkx_labels(ego, pos, labels=labels, font_size=12, font_weight='bold')
        nx.draw_networkx_edges(ego, pos, width=[w / 5 for w in width.values()], edge_color=edge_colors)
#        nx.draw_networkx_edge_labels(ego, pos, node_size=sizes, node_color=color_map, edge_labels=width)
        plt.axis('off')
        plt.tight_layout()
        if save_fig: plt.savefig(file_name + ".pdf")
        plt.close()
        return ego, labels, pos, sizes, width, color_map, node_borders
    else:
        print(central_node, 'not in the network')
        return None, None, None, None, None, None, None
    
def create_perturbation_network(net1, net2, db1, db2, threshold=4):
    # check which edges are only in net1 or only in net2 and which are in common:
    edges1 = set(net1.edges)
    edges2 = set(net2.edges)
    common_edges = edges1.intersection(edges2)
    edges_1only = edges1.difference(edges2)
    edges_2only = edges2.difference(edges1)
    
    net = nx.Graph()
    edge_colors = {}
    size_map = []
    node_borders = []
    color_map = []
    node_labels = {}
    added_edges = []
    
    for u, v in common_edges:
        wij1 = net1.get_edge_data(u, v)['weight']
        wij2 = net2.get_edge_data(u, v)['weight']
        deltawij = wij2 - wij1
        if deltawij > threshold:
            color = 'green'
            net.add_edge(u, v, weight=abs(deltawij))
            edge_colors[(u, v)] = color
            added_edges.append((u, v))
        elif deltawij < -threshold:
            color = 'red'
            net.add_edge(u, v, weight=abs(deltawij))
            edge_colors[(u, v)] = color
            added_edges.append((u, v))
            
    for u, v in edges_1only:
        wij = net1.get_edge_data(u, v)['weight']
        color = 'red'
        if wij > threshold:
            net.add_edge(u, v, weight=wij)
            edge_colors[(u, v)] = color
            added_edges.append((u, v))
        
    for u, v in edges_2only:
        wij = net2.get_edge_data(u, v)['weight']
        color = 'green'
        if wij > threshold:
            net.add_edge(u, v, weight=wij)
            edge_colors[(u, v)] = color
            added_edges.append((u, v))
    
    for node in net.nodes:
        k1 = net1.degree(node)
        k2 = net2.degree(node)
        w1 = net1.degree(node, weight='weight')
        w2 = net2.degree(node, weight='weight')
        if isinstance(k1, int) and isinstance(k2, int):
            nw1 = w1 / k1
            nw2 = w2 / k2
            
            info1 = db1[db1['Position'] == node]
            if len(info1) > 1: info1 = info1.iloc[0] # check why more than one
            type_aa_1 = amino_acids_conversion.three2one(info1["Type of residue"].item())
            
            info2 = db2[db2['Position'] == node]
            if len(info2) > 1: info2 = info2.iloc[0] # check why more than one
            type_aa_2 = amino_acids_conversion.three2one(info2["Type of residue"].item())
            
            if type_aa_1 == type_aa_2:
                color_map.append('gray')
                label = type_aa_1 + node[1::] + ":" + node[0]
                node_labels[node] = label
            
            else:
                color_map.append('blue')
                label = type_aa_1 + '-' + type_aa_2 + node[1::] + ":" + node[0]
                node_labels[node] = label
                
        elif isinstance(k2, int):
            nw1 = 0
            nw2 = w2 / k2
            color_map.append('green')
            info2 = db2[db2['Position'] == node]
            if len(info2) > 1: info2 = info2.iloc[0] # check why more than one
            type_aa_2 = amino_acids_conversion.three2one(info2["Type of residue"].item())
            label = type_aa_2 + node[1::] + ":" + node[0]
            node_labels[node] = label
            
        else:
            nw1 = w1 / k1
            nw2 = 0
            color_map.append('red')
            info1 = db1[db1['Position'] == node]
            if len(info1) > 1: info1 = info1.iloc[0] # check why more than one
            type_aa_1 = amino_acids_conversion.three2one(info1["Type of residue"].item())
            label = type_aa_1 + node[1::] + ":" + node[0]
            node_labels[node] = label
            
        deltanw = nw2 - nw1
        size_map.append(abs(deltanw)* 100 + 500)
        if deltanw > 0:
            node_borders.append('green')
        elif deltanw < 0:
            node_borders.append('red')
        else:
            node_borders.append('gray')
#    
#    node_labels = {node: node for node in net.nodes} #TO DO: labels as in normal networks, but with mutations
    
#    edge_color_map = [edge_colors[edge] for edge in added_edges]

    return net, node_labels, size_map, color_map, edge_colors, node_borders, edges_1only, edges_2only, common_edges

def create_ego_perturbation_network(net, central_node, all_labels, size_map,
                                    color_map, node_borders, edges_color_map):
        ego = nx.ego_graph(net, central_node)
        labels = {node: all_labels[node] for node in ego.nodes}
        neighbors = ego.copy()
        neighbors.remove_node(central_node)
        
        node_borders_dict = {node: node_borders[i] for i, node in enumerate(net.nodes)}
        color_map_dict = {node: color_map[i] for i, node in enumerate(net.nodes)}
        size_map_dict = {node: size_map[i] for i, node in enumerate(net.nodes)}
        
        node_borders_ego = {node: node_borders_dict[node] for node in ego.nodes}
        color_map_ego = [color_map_dict[node] for node in ego.nodes]
        size_map_ego = [size_map_dict[node] for node in ego.nodes]
        
        edges_color_map_dict = {edge: edges_color_map[i] for i, edge in enumerate(net.edges)}
        edges_color_map_ego = []
        for u, v in ego.edges:
            try:
                edges_color_map_ego.append(edges_color_map_dict[(u, v)])
            except KeyError:
                edges_color_map_ego.append(edges_color_map_dict[(v, u)])
#        edges_color_map_ego = {edge: edges_color_map_dict[edge] for edge in ego.edges}
        
        if not node_borders:
            neighbors_structure = {}
            node_borders = {}
            node_borders[central_node] = 'grey'
            for neighbor in neighbors:
                struct, dist = type_structure(central_node, neighbor) # to modify
                neighbors_structure[neighbor] = (struct, dist)
                if struct == '1D': node_borders[neighbor] = 'turquoise'
                elif struct == '2D': node_borders[neighbor] = 'turquoise'
                elif struct == '3D': node_borders[neighbor] = 'cadetblue'
                else: node_borders[neighbor] = 'cadetblue'
        
#        else:
#            neighbors_structure = {neighbor: (type_structure(central_node, neighbor)) for neighbor in neighbors}
#            print(neighbors_structure)
#                    
#        sorted_neighbors = deque(sorted(neighbors_structure,
#                                        key=lambda x: neighbors_structure[x][1]))
#        print(sorted_neighbors)
#        initial = central_node[0] + str(int(central_node[1::]) + 1)
#        if initial in sorted_neighbors:
#            sorted_neighbors.rotate(-sorted_neighbors.index(initial))
#        else:
#            initial = central_node[0] + str(int(central_node[1::]) - 1)
#            sorted_neighbors.rotate(-sorted_neighbors.index(initial) + 1)
#        
#        pos_original=nx.circular_layout(neighbors)
#        pos = {}
#        for i in range(len(sorted_neighbors)):
#            pos[sorted_neighbors[i]] = list(pos_original.values())[i]
#        pos[central_node] = np.array([0, 0])
#        width = nx.get_edge_attributes(ego, 'weight')

        pos_original=nx.circular_layout(neighbors)
        pos = {}
        for i in range(len(neighbors.nodes)):
            pos[list(neighbors.nodes)[i]] = list(pos_original.values())[i]
        pos[central_node] = np.array([0, 0])
        width = nx.get_edge_attributes(ego, 'weight')
        
        return ego, labels, pos, size_map_ego, width, color_map_ego, node_borders_ego, edges_color_map_ego

def draw_neighborhood_perturbation(net, pos, all_labels, size_map, color_map, edge_color_map, node_borders, save_fig=True,
                      file_name='figure'):
    plt.figure(figsize=(6, 6))
    try:
        int(pos[0])
        central_node = 'A' + pos
    except ValueError:
        central_node = pos
    if central_node in net.nodes:
        ego, labels, pos, sizes, width, color_map, node_borders, edges_color_map_ego = create_ego_perturbation_network(net, central_node, all_labels, size_map,
                                    color_map, node_borders, edge_color_map)
        nx.draw_networkx_nodes(ego, pos, node_size=sizes, node_color=color_map, edgecolors=[node_borders[n] for n in ego.nodes], linewidths=4)
        nx.draw_networkx_labels(ego, pos, labels=labels, font_size=12, font_weight='bold')
        nx.draw_networkx_edges(ego, pos, width=[w / 5 for w in width.values()], edge_color=edges_color_map_ego)
#        nx.draw_networkx_edge_labels(ego, pos, node_size=sizes, node_color=color_map, edge_labels=width)
        plt.axis('off')
        plt.tight_layout()
        if save_fig: plt.savefig(file_name + ".pdf")
        plt.close()
        return ego, labels, pos, sizes, width, color_map, node_borders, edge_color_map
    else:
        print(central_node, 'not in the network')
        return None, None, None, None, None, None, None