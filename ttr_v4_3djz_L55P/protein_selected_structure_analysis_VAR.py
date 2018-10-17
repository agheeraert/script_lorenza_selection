import sys
import os
from glob import glob
from Bio.PDB import PDBList
import networkx as nx
sys.path.append(os.path.realpath(__file__).rsplit('\\', 2)[0])
import biographs as bg
from secondary_structure import assign_secondary_structure, get_neighbor_structure_relation
import csv
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from networkx.algorithms import community
import numpy as np
import network_visualization
import pickle

selected_positions = range(11, 125)
reference_folder = 'ttr_v4_1f41_WT_HUM'
dim = '3-4D'
folder_path = dim + 'analysis'
source = 'A55'

if dim == '2D':
    rel_list = ['2D2', '2D3', '2D4']
elif dim == '1D':
    rel_list = ['2D1']
elif dim == '3-4D':
    rel_list = ['3D', '4D']
else:
    rel_list = [dim]

plt.rc('text', usetex=True)
os.makedirs(folder_path, exist_ok=True)

if "single_mutations.csv" in os.listdir() and "single_mutations_non_pathogenic.csv" in os.listdir():
    mutation_analysis = True
else:
    mutation_analysis = False
    pathogenic = []
    non_pathogenic = []
    both = []
    
os.makedirs(folder_path + '\\neighborhoods', exist_ok=True)

# load k_w dictionary (from database)
cwd = os.getcwd()
db_dir = cwd.rsplit('\\', 1)[0]
k_w = pickle.load(open(db_dir + '\\k_w.p', 'rb'))

# CREATE DATABASE

current_path = os.path.realpath(__file__).rsplit('\\', 1)[0]

with open('aminoacids.txt', 'r') as f:
    amino_acids = [aa.rsplit('\n')[0] for aa in f]

# Create a temporal directory called `pdb'
pdbs_path = os.path.join(current_path, 'pdb')

if pdbs_path not in glob(os.path.join(current_path, '*')):
    os.mkdir(pdbs_path)

with open(os.path.join(current_path, 'pdbs.txt'), 'r') as f:
    pdbs = [pdb[:-1] for pdb in f]
    if len(pdbs) == 1:
        pdb_id = pdbs[0]
    else:
        pdb_id = current_path.rsplit('\\', 1)[1]

pdbl = PDBList(obsolete_pdb=True)

if not glob(os.path.join(pdbs_path, '*')):
    pdbl.download_pdb_files(pdbs, file_format='pdb', pdir=pdbs_path)

#initialize databases to report
database_1 = []
database_2 = []


pdbs = glob(os.path.join(pdbs_path, '*'))

for pdb in pdbs:
    mol = bg.Pmolecule(pdb)
    net = mol.network()
    
    # take only selected positions:
    if selected_positions:
        for node in list(net.nodes):
            pos = int(node[1::])
            if pos not in selected_positions:
                net.remove_node(node)
    
    secondary_structure = assign_secondary_structure(pdb)
    
    residues_dict = {}
    for residue in mol.model.get_residues():
        res_type = residue.resname
        res_pos = residue.parent.id + str(residue.id[1])
        residues_dict[res_pos] = res_type
    
    for residue in mol.model.get_residues():
        adj_vector = [0] * 20
        weight_vector = [0] * 20
        node_name = residue.parent.id + str(residue.id[1])
        deg = nx.degree(net, residue.parent.id + str(residue.id[1]))
        if deg == 0:
            net.remove_node(residue.parent.id + str(residue.id[1]))
        else:
            weight = nx.degree(net, residue.parent.id + str(residue.id[1]),
                               weight='weight')
            restype = residue.resname
            resname = pdb.rsplit('\\', 1)[1][:-4] + residue.parent.id \
                    + str(residue.id[1])
            size = len(residue)
            seqpos = residue.id[1]
            if seqpos not in selected_positions:
                continue
            structure = secondary_structure[node_name]
            
            # check how many other aas can have the same k and w in the database
            n_others = 0
            for aa in amino_acids:
                if aa != restype:
                    try:
                        [w_min, w_max] = k_w[aa][deg]
                        if w_min <= weight and w_max >= weight:
                            n_others += 1
                    except KeyError:
                        pass
                    
            line = [resname, node_name, seqpos, restype, deg, weight, weight/deg,
                size, structure, n_others]
            line_2 = [resname, node_name, seqpos, restype, deg, weight, weight/deg,
                  size, structure, n_others]
            
            for neighbor in list(nx.neighbors(net, node_name)):
                neighbor_type = residues_dict[neighbor]
                edge_weight = nx.edges(net)[(node_name, neighbor)]['weight']
                aa_num = amino_acids.index(neighbor_type)
                adj_vector[aa_num] += 1
                weight_vector[aa_num] += edge_weight
                relation = get_neighbor_structure_relation(secondary_structure, node_name, neighbor)
                #select only 2D or 3D or 4D etc
                if relation in rel_list:
                    edge_name = resname[3:7] + '-' + node_name + '-' + neighbor
                    line.append(neighbor)
                    line.append(neighbor_type)
                    line.append(edge_weight)
                    line.append(relation)
                else:
                    net.remove_edge(neighbor, node_name)
            # check if the residue became of degree zero:
            deg = nx.degree(net, residue.parent.id + str(residue.id[1]))
            if deg == 0:
                net.remove_node(residue.parent.id + str(residue.id[1]))
            else:
                database_1.append(line)
                line_2 += adj_vector
                database_2.append(line_2)


sortedlist_pos = sorted(database_1, key=lambda row: row[2])

sortedlist_pos_2 = sorted(database_2, key=lambda row: row[2])

    
with open(folder_path + "\\database_pos_1.csv", "w", newline = '') as f:
    writer = csv.writer(f)
    writer.writerow(['Residue name','Position', 'Sequence position', 'Type of residue',
                     'Degree', 'Weight', 'Weight/Degree', 'Atomic number', 'Secondary structure', 'N. others',
                     'Neighbor position', 'Neighbor type', 'Pairwise weight', 'Relation'])
    writer.writerows(sortedlist_pos)
    
with open(folder_path + "\\database_pos_2.csv", "w", newline = '') as f:
    writer = csv.writer(f)
    writer.writerow(['Residue name','Position', 'Sequence position', 'Type of residue',
                     'Degree', 'Weight', 'Weight/Degree',
                     'Atomic number', 'Secondary structure', 'N. others'] + amino_acids)
    writer.writerows(sortedlist_pos_2)

# database 1 has to have all rows of the same lenght to be read as a dataframe
lengths = [len(row) for row in database_1]
max_length = max(lengths)

db_1 = [['Residue name','Position', 'Sequence position', 'Type of residue',
                     'Degree', 'Weight', 'Weight/Degree', 'Atomic number', 'Secondary structure', 'N. others',
                     'Neighbor position', 'Neighbor type', 'Pairwise weight', 'Relation']]
missing_header = max_length - len(db_1[0])
for i in range(int(missing_header / 4)):
    db_1[0].append('Neighbor position')
    db_1[0].append('Neighbor type')
    db_1[0].append('Pairwise weight')
    db_1[0].append('Relation')
   
for row in database_1:
    missing = max_length - len(row)
    for i in range(int(missing)):
        row.append('-')
    db_1.append(row)

db_1 = pd.DataFrame(db_1[1::], columns=db_1[0])

db_2 = pd.DataFrame(database_2, columns=['Residue name','Position', 'Sequence position', 'Type of residue',
                     'Degree', 'Weight', 'Weight/Degree',
                     'Atomic number', 'Secondary structure', 'N. others'] + amino_acids)

# NEIGHBORHOOD WATCH DISTRIBUTION
    
neigh_watch = db_2['Weight/Degree'].values

max_value = max(neigh_watch)

def herfindhal_index(x):
    num = sum([v ** 2 for v in x])
    den = (sum(x)) ** 2
    H = num / den
    return H

plt.figure(figsize=(10,10))
plt.hist(neigh_watch, bins=range(0, int(max_value + 1) + 1), normed=True, rwidth=0.9)
plt.title('Average weight distribution', fontsize=22)
plt.xlabel('w/k', fontsize=22)
plt.ylabel('P(w/k)', fontsize=22)
H = herfindhal_index(neigh_watch)
N = len(neigh_watch)
t = '$\\newline$'.join(['', 'Herfindhal index:', 'H = %.3f' %(H), '1/H = %.1f' %(1 / H),
                 'N = %s' %(N), 'N. outliers = %.3f * N' %(1 - 1 / (H * N))])
plt.text(0.1, 0.15, t, verticalalignment='center',
        horizontalalignment='left', fontsize=18)
plt.savefig(folder_path + '\\neighborhood_watch_distribution1.pdf')


fig, ax1 = plt.subplots(figsize=(10,10))
plt.hist(neigh_watch, bins=range(0, int(max_value + 1) + 1), rwidth=0.8)
plt.title('Average node weight distribution', fontsize=22)
plt.ylabel('N(w/k)', fontsize=22)
plt.xlabel('w/k', fontsize=22)
plt.xticks(fontsize=22)
#plt.xlim(-1, int(max_value + 1) + 1)
left, bottom, width, height = [0.15, 0.65, 0.35, 0.1]
ax2 = fig.add_axes([left, bottom, width, height])
ax2.boxplot(neigh_watch, vert=False)
plt.yticks([])
plt.xticks(fontsize=18)
plt.xlabel('w/k', fontsize=18)
plt.xlim(-1, int(max_value + 1) + 1)
plt.savefig(folder_path + '\\neighborhood_watch_distribution2.pdf')
plt.show()

        
# DRAW PROTEIN NETWORK
        
os.makedirs(folder_path + '\\network_pictures', exist_ok=True)

(net,
node_labels,
sizes,
color_map, _) = network_visualization.create_network(pdb_id,
            database=db_2, net=net,
            colors='neighborhood_watch_sharp', sizes='degree')
(net2,
 node_labels2,
 sizes2,
 color_map2, _ )= network_visualization.create_network(pdb_id,
               database=db_2,net=net, colors=None, sizes=None)
(net3,
 node_labels3,
 sizes3, color_map3,
 edge_color_map) = network_visualization.create_network(pdb_id,
               database=db_2, net=net, colors='pairwise_sharp',
               sizes='neighborhood_watch_sharp')

sharp_sizes = []
for s in sizes:
    if s <= 5:
        sharp_sizes.append(1)
    elif s <= 13:
        sharp_sizes.append(2)
    else:
        sharp_sizes.append(3)

sizes = [s * 100 for s in sharp_sizes]

pos_spring = nx.spring_layout(net)
pos_circular = nx.circular_layout(net)


network_visualization.draw_network(net, node_labels, color_map=color_map,
                                   draw_edges=False, draw_labels=False,
                                   name=folder_path + '\\network_pictures\\' + pdb_id + '_color', 
                                   figsize=(50, 50), pos=pos_spring)
network_visualization.draw_network(net, node_labels, sizes=sizes,
                                   draw_labels=False,
                                   name=folder_path + '\\network_pictures\\' + pdb_id + '_size',
                                   figsize=(50, 50), pos=pos_spring)
network_visualization.draw_network(net, node_labels, sizes=sizes,
                                   color_map=color_map, draw_labels=False,
                                   name=folder_path + '\\network_pictures\\' + pdb_id + '_colorsize',
                                   figsize=(50, 50), pos=pos_spring)

network_visualization.draw_network(net2, node_labels2, color_map='deepskyblue',
                                   draw_edges=False, draw_labels=False,
                                   name=folder_path + '\\network_pictures\\' + pdb_id + '_nodes',
                                   figsize=(50, 50), pos=pos_spring)
network_visualization.draw_network(net2, node_labels2, sizes=sizes2,
                                   color_map='deepskyblue', draw_labels=False,
                                   name=folder_path + '\\network_pictures\\' + pdb_id + '_links',
                                   figsize=(50, 50), pos=pos_spring)
network_visualization.draw_network(net3, node_labels3, sizes=[s / 50 for s in sizes3],
                                   edge_color_map=edge_color_map,
                                   draw_labels=True,
                                   name=folder_path + '\\network_pictures\\' + pdb_id +
                                   '_links_color', figsize=(50, 50),
                                   pos=pos_spring, facecolor='silver')

# LEGENDS
G1 = nx.Graph()
colors1 = ['blue', 'cyan', 'greenyellow', 'yellow', 'orange', 'red']
for n in range(1, 7):
    G1.add_node(str(n))

pos1 = {}
for i, node in enumerate(G1.nodes):
    pos1[node] = np.array([0, i + 0.15])

G2 = nx.Graph()
sizes2 = [1, 2, 3, 4]
sizes2 = [s * 400 for s in sizes2]
for n in range(1, 5):
    G2.add_node(str(n))

pos2 = {}
for i, node in enumerate(G2.nodes):
    pos2[node] = np.array([0, 3 * i + 0.15])
    
G2b = nx.Graph()
sizes2b = [1, 2, 3, 4, 5, 6]
sizes2b = [s * 400 for s in sizes2]
for n in range(1, 7):
    G2b.add_node(str(n))

pos2b = {}
for i, node in enumerate(G2b.nodes):
    pos2b[node] = np.array([0, 3 * i + 0.15])

plt.figure() 
nx.draw_networkx_nodes(G1, pos1, node_color=colors1,edgecolors='k')
plt.axis('off')
plt.xlim(-0.2, 3)
plt.ylim(-1, 6)
plt.rc('text', usetex=True)
plt.rc('font', family='calibri')
plt.text(0.5, 0, 'W/k $<$ 5',fontsize=16)
plt.text(0.5, 1, '5 $leq$ W/k $<$ 10',fontsize=16)
plt.text(0.5, 2, '10 $leq$ W/k $<$ 15',fontsize=16)
plt.text(0.5, 3, '15 $leq$ W/k $<$ 20',fontsize=16)
plt.text(0.5, 4, '20 $leq$ W/k $<$ 25',fontsize=16)
plt.text(0.5, 5, 'W/k $geq$ 25',fontsize=16)
plt.tight_layout()
plt.savefig(folder_path + '\\network_pictureslegend_color.pdf')
plt.close()

texts = ['$w_{ij}   <$ 10', '10 $leq  w_{ij}   <$ 20', '20 $leq  w_{ij}   <$ 30',
         '30 $leq  w_{ij}   <$ 40', '40 $leq  w_{ij}   <$ 50', '$w_{ij}   geq$ 50']
fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='calibri')
for index, c in enumerate(colors1):
    subf = int('61%s' %(index + 1))
    ax1 = fig.add_subplot(subf)
    ax1.add_patch(patches.Rectangle((0.1, 0.1), 0.1, 0.5, facecolor=c, edgecolor='k',
                                    linewidth=0.1))
    ax1.axis('off')
    ax1.text(0.3, 0, texts[index], fontsize=16)
plt.axis('off')
plt.tight_layout()
plt.savefig(folder_path + '\\network_pictureslegend_link_color.pdf')
plt.close()


plt.figure()  
nx.draw_networkx_nodes(G2, pos2, node_size=sizes2, node_color='w', edgecolors='k')
plt.axis('off')
plt.xlim(-0.3, 3)
plt.ylim(-1, 6)
plt.rc('text', usetex=True)
plt.rc('font', family='calibri')
plt.text(0.5, 0, 'k $<$ 5',fontsize=16)
plt.text(0.5, 1, '5 $leq$ k $<$ 9',fontsize=16)
plt.text(0.5, 2, '9 $leq$ k $<$ 13',fontsize=16)
plt.text(0.5, 3, 'k $geq$ 13',fontsize=16)
plt.tight_layout()
plt.savefig(folder_path + '\\network_pictureslegend_size_k.pdf')
plt.close()

plt.figure()  
nx.draw_networkx_nodes(G2b, pos2b, node_size=sizes2b, node_color='w', edgecolors='k')
plt.axis('off')
plt.xlim(-0.3, 3)
plt.ylim(-1, 17)
plt.rc('text', usetex=True)
plt.rc('font', family='calibri')
plt.text(0.5, 0, 'W/k $<$ 5',fontsize=16)
plt.text(0.5, 3, '5 $leq$ W/k $<$ 10',fontsize=16)
plt.text(0.5, 6, '10 $leq$ W/k $<$ 15',fontsize=16)
plt.text(0.5, 9, '15 $leq$ W/k $<$ 20',fontsize=16)
plt.text(0.5, 12, '20 $leq$ W/k $<$ 25',fontsize=16)
plt.text(0.5, 15, 'W/k $geq$ 25',fontsize=16)
plt.tight_layout()
plt.savefig(folder_path + '\\network_pictureslegend_size_nw.pdf')
plt.close()

# DRAW NEIGHBORHOODS

os.makedirs(folder_path + '\\neighborhoods', exist_ok=True)

for node in net.nodes:
    if node[0] != 'A': break            
    network_visualization.draw_neighborhood(net, node, node_labels, pathogenic,
                                            non_pathogenic, both, save_fig=True,
                                            file_name=folder_path + '\\neighborhoods\\'+
                                            pdb_id + node, threshold_edge=20,
                                            sort=False)


# DRAW PERTURBATION NETWORK
path_ref = os.getcwd().rsplit('\\', 1)[0] + '\\' + reference_folder
net_ref = pickle.load(open(path_ref + '\\' + dim + 'net_linkcolor.p', 'rb'))
db_ref = pd.DataFrame(pd.read_csv(path_ref + "\\database_pos_2.csv"))

(net_p, node_labels_p, size_map_p,
 color_map_p, edge_color_dict_p,
 node_borders_p, edges_1only,
 edges_2only,
 common_edges) = network_visualization.create_perturbation_network(net_ref,
               net3, db_ref, db_2)

edge_color_map_p = []
for u, v in net_p.edges:
    try:
        edge_color_map_p.append(edge_color_dict_p[u, v])
    except KeyError:
        edge_color_map_p.append(edge_color_dict_p[v, u])

if dim == '4D' or dim == '3-4D':
    nodes_pos = {}
    position_A = 0
    position_B = 0
    shift_A = 1
    shift_B = 1
    for i, node in enumerate(sorted(net_p.nodes, key= lambda x: int(x[1::]))):
        chain = node[0]
        if chain == 'A':
            nodes_pos[node] = (shift_A, position_A)
            position_A += 1
            shift_A *= -1
        else:
            nodes_pos[node] = (10 + shift_B, position_B)
            position_B += 1
            shift_B *= -1

else:
    nodes_pos = None
        
network_visualization.draw_network(net_p, node_labels_p, sizes=size_map_p,
                                   edge_color_map=edge_color_map_p,
                                   color_map=color_map_p,
                                   draw_labels=True,
                                   name= folder_path + '\\network_pictures\\' +
                                   pdb_id + 'perturbation_net_links_color',
                                   figsize=(20, 30),
                                   node_edge_colors=node_borders_p,
                                   pos = nodes_pos)

os.makedirs(folder_path + '\\neighborhoods_perturbation', exist_ok=True)

for node in net_p.nodes:
    if node[0] != 'A': continue
    network_visualization.draw_neighborhood_perturbation(net_p,
                                node, node_labels_p, size_map_p, color_map_p,
                                edge_color_map_p, node_borders_p,
                                file_name=folder_path + '\\neighborhoods_perturbation\\'+
                                pdb_id + node)

sizes_dict = {node: size_map_p[i] for i, node in enumerate(net_p.nodes)}
edge_color_dict = {edge: edge_color_map_p[i] for i, edge in enumerate(net_p.edges)}
color_dict = {node: color_map_p[i] for i, node in enumerate(net_p.nodes)}
node_borders_dict = {node: node_borders_p[i] for i, node in enumerate(net_p.nodes)}
weights_p = nx.get_edge_attributes(net_p, 'weight')

# PERTURBATION TREE

tree = nx.bfs_tree(net_p, source)

node_labels_tree = {node: node_labels_p[node] for node in tree.nodes}

sizes_tree = [sizes_dict[node] for node in tree.nodes]

#add triangles
for v in tree.nodes:
    for n in net_p.neighbors(v):
        if n in tree.nodes:
            tree.add_edge(v, n)

tree = tree.to_undirected()

edge_color_tree = []
weights_tree = {}
for u, v in tree.edges:
    try:
        edge_color_tree.append(edge_color_dict[(u, v)])
    except KeyError:
        edge_color_tree.append(edge_color_dict[(v, u)])
    try:
        weights_tree[(u, v)] = weights_p[(u, v)]
    except KeyError:
        weights_tree[(u, v)] = weights_p[(v, u)]

nx.set_edge_attributes(tree, weights_tree, name='weight')

color_tree = [color_dict[node] for node in tree.nodes]
node_borders_tree = [node_borders_dict[node] for node in tree.nodes]

pickle.dump([tree, node_labels_tree, sizes, edge_color_tree,
             color_tree, node_borders_tree], open('tree.p', 'wb'))
 
if dim == '4D' or dim == '3-4D':
    pos_tree = {}
    position_A = 0
    position_B = 0
    shift_A = 1
    shift_B = 1
    for i, node in enumerate(sorted(tree.nodes, key= lambda x: int(x[1::]))):
        chain = node[0]
        if chain == 'A':
            pos_tree[node] = (shift_A, position_A)
            position_A += 1
            shift_A *= -1
        else:
            pos_tree[node] = (3 + shift_B, position_B)
            position_B += 1
            shift_B *= -1

else:
    pos_tree = None
        
node_borders_tree = ['r']
node_borders_tree = node_borders_tree + ['gray' for node in list(tree.nodes)[1::]]

color_tree = ['y']
color_tree = color_tree + ['lightgray' for node in list(tree.nodes)[1::]]


plt.rc('text', usetex=False)
network_visualization.draw_network(tree, node_labels_tree, sizes=600,
                                   edge_color_map=edge_color_tree,
                                   color_map = color_tree,
                                   draw_labels=True,
                                   name= folder_path + '\\network_pictures\\' + pdb_id +
                                   'perturbation_tree_' + source + '_' + dim, figsize=(10, 10),
                                   pos=pos_tree, labels_size=18,
                                   node_edge_size=1)
