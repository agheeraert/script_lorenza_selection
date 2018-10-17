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
import numpy as np
import network_visualization
import pickle

selected_positions = range(11, 125) # consider only positions that are
                                    # common in all structures
dim = '3-4D'
folder_path = dim + 'analysis'

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
    
os.makedirs(folder_path + '\\neighborhoods', exist_ok=True)

if "single_mutations.csv" in os.listdir() and "single_mutations_non_pathogenic.csv" in os.listdir():
    mutation_analysis = True
else:
    mutation_analysis = False
    pathogenic = []
    non_pathogenic = []
    both = []
    
os.makedirs('neighborhoods', exist_ok=True)

# load k_w dictionary (from database)
cwd = os.getcwd()
db_dir = cwd.rsplit('\\', 1)[0]
k_w = pickle.load(open(db_dir + '\\k_w.p', 'rb'))

# CREATE AA NETWORK

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

# NOTE: I CALL IT DATABASE BUT IT CONCERNS ONLY THIS STRUCTURE.
# IT IS JUST DUE TO COPY-PASTE FROM THE DATABASE ANALYSIS. TO CORRECT.
    
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
                #select only 2D
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

columns = ['Residue name','Position', 'Sequence position', 'Type of residue',
             'Degree', 'Weight', 'Weight/Degree', 'Atomic number',
             'Secondary structure', 'N. others', 'Neighbor position',
             'Neighbor type', 'Pairwise weight', 'Relation']

db_1 = []
missing_header = max_length - len(columns)
for i in range(int(missing_header / 4)):
    columns.append('Neighbor position')
    columns.append('Neighbor type')
    columns.append('Pairwise weight')
    columns.append('Relation')
   
for row in database_1:
    missing = max_length - len(row)
    for i in range(int(missing)):
        row.append('-')
    db_1.append(row)

db_1 = pd.DataFrame(db_1, columns=columns)

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

nw_helix = []
nw_sheet = []
nw_loop = []

for _, row in db_2.iterrows():
    nw = row['Weight/Degree']
    secstr = row['Secondary structure']
    if secstr[0] == 'h':
        nw_helix.append(nw)
    elif secstr[0] == 's':
        nw_sheet.append(nw)
    else:
        nw_loop.append(nw)


plt.figure(figsize=(10,10))   
plt.hist([nw_helix, nw_sheet, nw_loop], bins=range(0, int(max_value + 1) + 1),
         label=['helix', 'sheet', 'loop'], alpha=1, histtype='barstacked',
         rwidth=0.9)
plt.legend()
plt.title("Node's average link weight distribution - database")
plt.xlabel('w/k')
plt.ylabel('P(w/k)')
H = herfindhal_index(neigh_watch)
N = len(neigh_watch)
t = '$\\newline$'.join(['', 'Herfindhal index:', 'H = %.3f' %(H), '1/H = %.1f' %(1 / H),
                 'N = %s' %(N), 'N. outliers = %.3f * N' %(1 - 1 / (H * N))])
plt.text(0.1, 25, t, verticalalignment='center',
        horizontalalignment='left', fontsize=18)
plt.savefig(folder_path + '\\neighborhood_watch_distribution1.pdf')


fig, ax1 = plt.subplots(figsize=(10,10))
plt.hist([nw_helix, nw_sheet, nw_loop], bins=range(0, int(max_value + 1) + 1),
         label=['helix', 'sheet', 'loop'], alpha=1, histtype='barstacked',
         rwidth=0.9)
plt.legend()
plt.title("Node's average link weight distribution - database")
plt.xlabel('w/k')
plt.ylabel('P(w/k)')
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

net, node_labels, sizes, color_map, _ = network_visualization.create_network(
        pdb_id, database=db_2, net=net, colors='neighborhood_watch_sharp',
        sizes='degree')
net2, node_labels2, sizes2, color_map2, _ = network_visualization.create_network(
        pdb_id, database=db_2, net=net, colors=None,
        sizes=None)
net3, node_labels3, sizes3, color_map3, edge_color_map = network_visualization.create_network(
        pdb_id, database=db_2, net=net, colors='pairwise_sharp',
        sizes='neighborhood_watch_sharp')

sharp_sizes = []
for s in sizes:
    if s <= 5:
        sharp_sizes.append(1)
    elif s <= 13:
        sharp_sizes.append(2)
    else:
        sharp_sizes.append(3)

sizes = [s * 2000 for s in sharp_sizes]

pos_spring = nx.spring_layout(net)
pos_circular = nx.circular_layout(net)


network_visualization.draw_network(net,
                                   node_labels, color_map=color_map,
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
network_visualization.draw_network(net3, node_labels3, sizes=sizes3,
                                   edge_color_map=edge_color_map,
                                   draw_labels=False,
                                   name=folder_path + '\\network_pictures\\' + pdb_id +
                                   '_links_color', figsize=(50, 50),
                                   pos=pos_spring, facecolor='silver')


pickle.dump(net, open(dim + 'net_nodecolor.p', 'wb'))
pickle.dump(net2, open(dim + 'net_nocolor.p', 'wb'))
pickle.dump(net3, open(dim + 'net_linkcolor.p', 'wb'))

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
plt.text(0.5, 1, '5 $\leq$ W/k $<$ 10',fontsize=16)
plt.text(0.5, 2, '10 $\leq$ W/k $<$ 15',fontsize=16)
plt.text(0.5, 3, '15 $\leq$ W/k $<$ 20',fontsize=16)
plt.text(0.5, 4, '20 $\leq$ W/k $<$ 25',fontsize=16)
plt.text(0.5, 5, 'W/k $\geq$ 25',fontsize=16)
plt.tight_layout()
plt.savefig(folder_path + '\\network_pictures\legend_color.pdf')
plt.close()

texts = ['$w_{ij} \  <$ 10', '10 $\leq \ w_{ij} \  <$ 20', '20 $\leq \ w_{ij} \  <$ 30',
         '30 $\leq \ w_{ij} \  <$ 40', '40 $\leq \ w_{ij} \  <$ 50', '$w_{ij} \  \geq$ 50']
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
plt.savefig(folder_path + '\\network_pictures\legend_link_color.pdf')
plt.close()


plt.figure()  
nx.draw_networkx_nodes(G2, pos2, node_size=sizes2, node_color='w', edgecolors='k')
plt.axis('off')
plt.xlim(-0.3, 3)
plt.ylim(-1, 6)
plt.rc('text', usetex=True)
plt.rc('font', family='calibri')
plt.text(0.5, 0, 'k $<$ 5',fontsize=16)
plt.text(0.5, 1, '5 $\leq$ k $<$ 9',fontsize=16)
plt.text(0.5, 2, '9 $\leq$ k $<$ 13',fontsize=16)
plt.text(0.5, 3, 'k $\geq$ 13',fontsize=16)
plt.tight_layout()
plt.savefig(folder_path + '\\network_pictures\legend_size_k.pdf')
plt.close()

plt.figure()  
nx.draw_networkx_nodes(G2b, pos2b, node_size=sizes2b, node_color='w', edgecolors='k')
plt.axis('off')
plt.xlim(-0.3, 3)
plt.ylim(-1, 17)
plt.rc('text', usetex=True)
plt.rc('font', family='calibri')
plt.text(0.5, 0, 'W/k $<$ 5',fontsize=16)
plt.text(0.5, 3, '5 $\leq$ W/k $<$ 10',fontsize=16)
plt.text(0.5, 6, '10 $\leq$ W/k $<$ 15',fontsize=16)
plt.text(0.5, 9, '15 $\leq$ W/k $<$ 20',fontsize=16)
plt.text(0.5, 12, '20 $\leq$ W/k $<$ 25',fontsize=16)
plt.text(0.5, 15, 'W/k $\geq$ 25',fontsize=16)
plt.tight_layout()
plt.savefig(folder_path + '\\network_pictures\legend_size_nw.pdf')
plt.close()


# DRAW NEIGHBORHOODS
for node in net.nodes:
    if node[0] != 'A': break            
    network_visualization.draw_neighborhood(net, node, node_labels, pathogenic,
                                            non_pathogenic, both, save_fig=True,
                                            file_name=folder_path + '\\neighborhoods\\'+ pdb_id + node,
                                            threshold_edge=20)
