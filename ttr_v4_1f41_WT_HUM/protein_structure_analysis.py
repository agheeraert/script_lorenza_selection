import sys
import os
from glob import glob
from Bio.PDB import PDBList
import networkx as nx
sys.path.append(os.path.realpath(__file__).rsplit('/', 2)[0])
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

plt.rc('text', usetex=True)

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
db_dir = cwd.rsplit('/', 1)[0]
k_w = pickle.load(open(db_dir + '/k_w.p', 'rb'))

# CREATE AA NETWORK

current_path = os.path.realpath(__file__).rsplit('/', 1)[0]

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
        pdb_id = current_path.rsplit('/', 1)[1]

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
        weight = nx.degree(net, residue.parent.id + str(residue.id[1]),
                           weight='weight')
        restype = residue.resname
        resname = pdb.rsplit('/', 1)[1][:-4] + residue.parent.id \
                + str(residue.id[1])
        size = len(residue)
        seqpos = residue.id[1]
        if seqpos not in selected_positions:
            continue
        structure = secondary_structure[node_name]
        
        # check how many others aas can have the same k and w in the database
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
            edge_name = resname[3:7] + '-' + node_name + '-' + neighbor
#            wijdict[restype][neighbor_type][relation] = wijdict[restype][neighbor_type][relation] + [(edge_name, edge_weight)]
            line.append(neighbor)
            line.append(neighbor_type)
            line.append(edge_weight)
            line.append(relation)
        database_1.append(line)
        line_2 += adj_vector
        database_2.append(line_2)


sortedlist_pos = sorted(database_1, key=lambda row: row[2])

sortedlist_pos_2 = sorted(database_2, key=lambda row: row[2])

    
with open("database_pos_1.csv", "w", newline = '') as f:
    writer = csv.writer(f)
    writer.writerow(['Residue name','Position', 'Sequence position', 'Type of residue',
                     'Degree', 'Weight', 'Weight/Degree', 'Atomic number', 'Secondary structure', 'N. others',
                     'Neighbor position', 'Neighbor type', 'Pairwise weight', 'Relation'])
    writer.writerows(sortedlist_pos)
    
with open("database_pos_2.csv", "w", newline = '') as f:
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
t = '$/newline$'.join(['', 'Herfindhal index:', 'H = %.3f' %(H), '1/H = %.1f' %(1 / H),
                 'N = %s' %(N), 'N. outliers = %.3f * N' %(1 - 1 / (H * N))])
plt.text(0.1, 25, t, verticalalignment='center',
        horizontalalignment='left', fontsize=18)
plt.savefig('neighborhood_watch_distribution1.pdf')


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
plt.savefig('neighborhood_watch_distribution2.pdf')
plt.show()

# NEIGHBORHOOD WATCH ALONG THE SEQUENCE

nw_chain = []
seq_positions = []
# select only chain A
for _, row in db_2.iterrows():
    if row['Position'][0] == 'B':
        break
    else:
        nw_chain.append(row['Weight/Degree'])
        seq_positions.append(row['Sequence position'])
        
colors = ['b' if nw >= 10 and nw <= 16 else 'r' for nw in nw_chain]

plt.figure(figsize=(35, 3.5))
plt.scatter(seq_positions, nw_chain, c=colors)
plt.ylabel('w/k')
plt.xlabel('sequence')
plt.xticks(seq_positions)
plt.yticks([min(nw_chain), 10, 16, max(nw_chain)])
plt.xlim(min(seq_positions) - 0.5, max(seq_positions) + 0.5)
plt.ylim(6, 25)
plt.grid()
plt.savefig('nw_sequence.pdf')
plt.show()

#if mutation_analysis:
#    # LOAD MUTATIONS
#    mutations_type = {}
#    mutations = []
#    pathogenic = []
#    pathogenic_names = []
#    with open("single_mutations.csv", "r") as f:
#        reader = csv.reader(f)
#        for line in reader:
#            pathogenic.append(line[1])
#            name = ''.join(line).replace(' ', '')
#            pathogenic_names.append(name)
#            mutations_type[name] = 'pathogenic'
#            mutations.append(line)
#    
#    non_pathogenic = []
#    non_pathogenic_names = []
#    with open("single_mutations_non_pathogenic.csv", "r") as f:
#        reader = csv.reader(f)
#        for line in reader:
#            non_pathogenic.append(line[1])
#            name = ''.join(line).replace(' ', '')
#            non_pathogenic_names.append(name)
#            mutations_type[name] = 'non_pathogenic'
#            mutations.append(line)
#    
#    def remove_duplicates(values):
#        output = []
#        seen = set()
#        for value in values:
#            # If value has not been encountered yet,
#            # ... add it to both list and set.
#            if value not in seen:
#                output.append(value)
#                seen.add(value)
#        return output
#    
#    pathogenic = remove_duplicates(pathogenic)
#    non_pathogenic = remove_duplicates(non_pathogenic)
#    
#    both = []
#    for p in pathogenic:
#        if p in non_pathogenic:
#            pathogenic.pop(pathogenic.index(p))
#            non_pathogenic.pop(non_pathogenic.index(p))
#            both.append(p)
#
#    # CHECK "WRONG" DEGREE, WEIGHT, NEIGHBORHOOD WATCH
#            
#columns = ['Position',
#           'Degree', 'Weight', 'Weight/Degree',
#           'Check Degree', 'Check Weight', 'Check Weight/Degree', 'Mutation?']
#    
#to_report = []
#
#for index, row in db_2.iterrows():
#    pos = row['Position']
#    k = row['Degree']
#    w = row['Weight']
#    nw = row['Weight/Degree']
#    mut =  ('Pathogenic' if pos[1::] in pathogenic else 'Non pathogenic' if pos[1::] in non_pathogenic else 'Both' if pos[1::] in both else '-')
#    info = [pos, k, w, nw,
#            True if k>= 4 and k <= 14 else False,
#            True if w>= 60 and w <= 140 else False,
#            True if nw>= 10 and nw <= 16 else False,
#            mut]
#    to_report.append(info)
#    
#report = pd.DataFrame(to_report, columns=columns)
#
#tot = len(report)
#
#wrong_k = tot - sum(report['Check Degree'])
#wrong_w = tot - sum(report['Check Weight'])
#wrong_nw = tot - sum(report['Check Weight/Degree'])
#
#print('wrong k: ', wrong_k, '/', tot)
#print('wrong w: ', wrong_w, '/', tot)
#print('wrong nw: ', wrong_nw, '/', tot)
#
#report.to_csv('check_deg_weight.csv', index=False)
#
#with open('check_deg_weight.csv', 'a', newline = '') as f:
#    writer = csv.writer(f)
#    writer.writerow(['N. WRONG', '-', '-', '-', str(wrong_k),
#                      str(wrong_w), str(wrong_nw), '-'])
        
# DRAW PROTEIN NETWORK

os.makedirs('network_pictures', exist_ok=True)

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
                                   name='network_pictures/' + pdb_id + '_color',
                                   figsize=(50, 50), pos=pos_spring)
network_visualization.draw_network(net, node_labels, sizes=sizes,
                                   draw_labels=False,
                                   name='network_pictures/' + pdb_id + '_size',
                                   figsize=(50, 50), pos=pos_spring)
network_visualization.draw_network(net, node_labels, sizes=sizes,
                                   color_map=color_map, draw_labels=False,
                                   name='network_pictures/' + pdb_id + '_colorsize',
                                   figsize=(50, 50), pos=pos_spring)

network_visualization.draw_network(net2, node_labels2, color_map='deepskyblue',
                                   draw_edges=False, draw_labels=False,
                                   name='network_pictures/' + pdb_id + '_nodes',
                                   figsize=(50, 50), pos=pos_spring)
network_visualization.draw_network(net2, node_labels2, sizes=sizes2,
                                   color_map='deepskyblue', draw_labels=False,
                                   name='network_pictures/' + pdb_id + '_links',
                                   figsize=(50, 50), pos=pos_spring)
network_visualization.draw_network(net3, node_labels3, sizes=sizes3,
                                   edge_color_map=edge_color_map,
                                   draw_labels=False,
                                   name='network_pictures/' + pdb_id +
                                   '_links_color', figsize=(50, 50),
                                   pos=pos_spring, facecolor='silver')


pickle.dump(net, open('net_nodecolor.p', 'wb'))
pickle.dump(net2, open('net_nocolor.p', 'wb'))
pickle.dump(net3, open('net_linkcolor.p', 'wb'))

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
plt.savefig('network_pictures\legend_color.pdf')
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
plt.savefig('network_pictures\legend_link_color.pdf')
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
plt.savefig('network_pictures\legend_size_k.pdf')
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
plt.savefig('network_pictures\legend_size_nw.pdf')
plt.close()


#if mutation_analysis:
#    # NETWORK MUTATIONS          
#            
#    net2, node_labels2, size_map2, color_map2, edge_color_map2 = network_visualization.create_network(pdb_id, database=db_2, net=net, colors='mutations', sizes='neighborhood_watch_sharp',
#                                                        pathogenic=pathogenic, non_pathogenic=non_pathogenic, both=both)
#    
#    network_visualization.draw_network(net, node_labels, sizes=sizes, color_map=color_map2, draw_labels=False, name='network_pictures/' + pdb_id + '_mutations', figsize=(50, 50), pos=pos_spring)
#
#
#    # LEGENDS
#    G3 = nx.Graph()
#    colors3 = ['gold', 'limegreen', 'tomato']
#    for n in range(1, 4):
#        G3.add_node(str(n))
#    
#    pos3 = {}
#    for i, node in enumerate(G3.nodes):
#        pos3[node] = np.array([0, i + 0.3])
#        
#    G4 = nx.Graph()
#    sizes4 = range(1, 7)
#    sizes4 = [s * 220 for s in sizes4]
#    for n in range(1, 7):
#        G4.add_node(str(n))
#    
#    pos4 = {}
#    for i, node in enumerate(G4.nodes):
#        pos4[node] = np.array([0, 4 * i + 0.5])
#    
#    plt.figure() 
#    nx.draw_networkx_nodes(G3, pos3, node_color=colors3,edgecolors='k')
#    plt.axis('off')
#    plt.xlim(-0.2, 3)
#    plt.ylim(-1, 3)
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='calibri')
#    plt.text(0.5, 0, 'Both',fontsize=16)
#    plt.text(0.5, 1, 'Non pathogenic',fontsize=16)
#    plt.text(0.5, 2, 'Pathogenic',fontsize=16)
#    plt.tight_layout()
#    plt.savefig('network_pictures\legend_color_mut.pdf')
#    plt.close()
#    
#    
#    plt.figure()  
#    nx.draw_networkx_nodes(G4, pos4, node_size=sizes4, node_color='w', edgecolors='k')
#    plt.axis('off')
#    plt.xlim(-0.5, 3)
#    plt.ylim(-2, 25)
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='calibri')
#    plt.text(0.5, 0, 'W/k $<$ 5',fontsize=16)
#    plt.text(0.5, 4, '5 $\leq$ W/k $<$ 10',fontsize=16)
#    plt.text(0.5, 8, '10 $\leq$ W/k $<$ 15',fontsize=16)
#    plt.text(0.5, 12, '15 $\leq$ W/k $<$ 20',fontsize=16)
#    plt.text(0.5, 16, '20 $\leq$ W/k $<$ 25',fontsize=16)
#    plt.text(0.5, 20, 'W/k $\leq$ 25',fontsize=16)
#    plt.tight_layout()
#    plt.savefig('network_pictures\legend_size_mut.pdf')
#    plt.close()  

    # NETWORK ANALYSIS
    
    # communities
    
    #communities_generator = community.girvan_newman(net)
    #top_level_communities = next(communities_generator)
    #next_level_communities = next(communities_generator)
    #third_level_communities = next(communities_generator)
    #
    #for i, comm in enumerate(third_level_communities):
    #    n_patho = 0
    #    n_non_patho = 0
    #    n_both = 0
    #    for node in comm:
    #        if node[1::] in pathogenic:
    #            n_patho += 1
    #        elif node[1::] in non_pathogenic:
    #            n_non_patho += 1
    #        elif node[1::] in both:
    #            n_both += 1
    #    print('-------')
    #    tot = len(comm)
    #    print('community ', i)
    #    print('patho: ', n_patho, n_patho/tot)
    #    print('non_patho: ', n_non_patho, n_non_patho/tot)
    #    print('both: ', n_both, n_both/tot)
    #    print('total: ', tot)
    
    # RESULT: patho e non patho homogeneously spread in all the communities
    
    
#    colors_dict = {}
#    n_colors = {}
#    for i, n in enumerate(net.nodes):
#        c = color_map2[i]
#        colors_dict[n] = c
#        try:
#            n_colors[c] = n_colors[c] + 1
#        except KeyError:
#            n_colors[c] = 1
#    
#    # number of links based on nodes' colors:
#    links_colors = {}
#    for (u, v) in net.edges:
#        cu = colors_dict[u]
#        cv = colors_dict[v]
#        try:
#            links_colors[(cu, cv)] = links_colors[(cu, cv)] + 1
#        except KeyError:
#            try:
#                links_colors[(cv, cu)] = links_colors[(cv, cu)] + 1
#            except KeyError:
#                links_colors[(cu, cv)] = 1
#                
#    # observed probability to peak a link connecting nodes of certain colors:
#    n_links = len(net.edges)
#    n_nodes = len(net.nodes)
#    
#    perc_link_colors = {key: value * 100 / n_links for key, value in links_colors.items()}
#    
#    prob_node_colors = {key: value / n_nodes for key, value in n_colors.items()}
#    
#    # hp: P(cu, cv) = P(cu)P(cv) -- note: I am not considering the degrees
#    theor_perc_link_colors = {}
#    for cu, cv in perc_link_colors.keys():
#        prob_cu = prob_node_colors[cu]
#        prob_cv = prob_node_colors[cv]
#        theor_perc_link_colors[(cu, cv)] = prob_cu * prob_cv * 100
#    
#    color_conversion = {
#            'tomato': 'pathogenical',
#            'lightgrey': 'nothing',
#            'limegreen': 'non pathogenical',
#            'gold': 'both'}
#    
#    to_report = [['color_1', 'color_2', 'observed probability', 'theoretical_probability']]
#    for cu, cv in perc_link_colors.keys():
#        to_report.append([
#                color_conversion[cu], color_conversion[cv],
#                perc_link_colors[(cu, cv)],
#                theor_perc_link_colors[(cu, cv)]])
#        
#    
#    with open('link_color_probability.csv', 'w', newline='') as f:
#        writer = csv.writer(f)
#        writer.writerows(to_report)
#     
#    nx.set_node_attributes(net, colors_dict, 'color')
#    a = nx.assortativity.attribute_assortativity_coefficient(net, 'color')
#    print('-------')
#    print('-------')
#    print('assortativity coefficient: ', a)
#    
#    # MUTATIONS ANALYSIS
#    
#    summary_patho = []
#    summary_non_patho = []
#    report_cases = {}
#    number_cases = {}
#    
#    wijdict = pickle.load(open(db_dir + '/wijdict.p', 'rb'))
#    
#    for [aa1, pos, aa2] in mutations:
#        name_mut = aa1 + pos + aa2
#        name_mut=name_mut.replace(' ', '')
#        report_cases[name_mut] = {}
#        number_cases[name_mut] = {}
#        check_aa1 = False
#        pos = pos.lstrip(' ')
#        
#        # take the information of the neighborhood in that position in the WT -- type2
#        row = db_2[db_2['Position'] == 'A' + pos].values
#        
#        # if the position in in the structure and the mutation is not a deletion:
#        if len(row) != 0 and aa2 != 'DEL':
#            # check that the 'initial' amino-acids is actually the one in the sequence at that position
#            if row[0][3] == aa1: check_aa1 = True
#            
#             # take the information of the neighborhood in that position in the WT -- type1
#            row1 = db_1[db_1['Position'] == 'A' + pos]
#            row1_name = row1['Residue name'].values[0]
#            k = int(row1['Degree'].values[0])
#            weight = int(row1['Weight'].values[0])
#            
#            # select the columns with info on the neighbors
#            neighbors_info = [list(row1.values[0][i:i+4]) for i in range(
#                    10, len(row1.values[0])-3, 4) if '-' not in list(row1.values[0][i:i+4])]
#        
#            # check in the database if every wij is possible after the mutation (tolerance of 6, check for the same 'structure relation')
#            for index, (neigh_pos, neigh_type, wij_wt, relation) in enumerate(neighbors_info):
#                cases = []
#                for name, wij_db in wijdict[aa2][neigh_type][relation]:
#                    diff_wij = np.abs(int(wij_db) - int(wij_wt))
#                    if diff_wij <= 6:
#                        cases.append((name, diff_wij))
#                report_cases[name_mut][neigh_pos] = cases
#                number_cases[name_mut][neigh_pos] = len(cases)
#    #        network_visualization.draw_neighborhood(net, pos, node_labels, pathogenic, non_pathogenic, both, save_fig=True, file_name='neighborhoods//' + name_mut, n_cases=number_cases[name_mut], threshold_edge=20)    
#    #        plt.close()
#                
#            mut_type = mutations_type[name_mut]
#            
#            if mut_type == 'pathogenic':
#                summary_patho.append([name_mut, check_aa1, k, weight, weight/k] + [
#                        item for sublist in [[key, value] for key, value in number_cases[name_mut].items()] for item in sublist])
#            else:
#                summary_non_patho.append([name_mut, check_aa1, k, weight, weight/k] + [
#                        item for sublist in [[key, value] for key, value in number_cases[name_mut].items()] for item in sublist])
#        
#        else:
#            mut_type = mutations_type[name_mut]
#            if mut_type == 'pathogenic':
#                summary_patho.append([name_mut] + ['-'] * 4)
#            else:
#                summary_non_patho.append([name_mut] + ['-'] * 4)
#                
#    with open("mutations_summary_patho.csv", 'w', newline = '') as f:
#        writer = csv.writer(f)
#        writer.writerow(['Mutation', 'Right initial aa?', 'Degree', 'Weight', 'Weight/Degree', 
#                         'Neighbor', 'N. cases same wij'])
#        writer.writerows(summary_patho)
#    
#    with open("mutations_summary_non_patho.csv", 'w', newline = '') as f:
#        writer = csv.writer(f)
#        writer.writerow(['Mutation', 'Right initial aa?', 'Degree', 'Weight', 'Weight/Degree', 
#                         'Neighbor', 'N. cases same wij'])
#        writer.writerows(summary_non_patho)
#    
#    
#    # CHECK "WRONG" PAIRWISE WEIGHTS
#    for aa in number_cases.keys():
#        if aa in pathogenic_names:
#            mut_type = 'patho'
#        else:
#            mut_type = 'non patho'
#        for aa2 in number_cases[aa].keys():
#            if number_cases[aa][aa2] < 20:
#                print(aa, aa2, number_cases[aa][aa2], mut_type)

# DRAW NEIGHBORHOODS
for node in net.nodes:
    if node[0] != 'A': break            
    network_visualization.draw_neighborhood(net, node, node_labels, pathogenic,
                                            non_pathogenic, both, save_fig=True,
                                            file_name='neighborhoods//'+ pdb_id + node,
                                            threshold_edge=20)

# Pairwise weight distributions

wij_list = []
for u, v in net.edges:
    wij = net.get_edge_data(u, v)['weight']
    wij_list.append(wij)

hist, bin_edges = np.histogram(wij_list, bins=range(1, max(wij_list) + 1), normed=True)
xticks = [x + 0.5 for x in bin_edges[0:-1]]

plt.figure()
plt.hist(wij_list, bins=range(1, max(wij_list) + 1, 3), normed=True)
#plt.scatter(xticks, hist)
plt.title('Pairwise weight distribution in the amino acids network')
plt.xlabel('w(i,j)')
plt.ylabel('P(w(i, j))')
#plt.ylim(min(hist), max(hist))
plt.yscale('log')
plt.savefig('pairwise_distrib.pdf')
plt.show()    


#subnetwork of w/k>13
subnet = net.copy()
lcc = True
to_remove = []
nw_subnet_size = []
for index, node in enumerate(subnet.nodes):
    k = nx.degree(subnet, node)
    w = nx.degree(subnet, node, weight='weight')
    nw = w / k
    if nw < 13:
        to_remove.append(node)
    else:
        if nw < 5: nw_subnet_size.append(1000)
        elif nw < 10: nw_subnet_size.append(4000)
        elif nw < 15: nw_subnet_size.append(7000)
        elif nw < 20: nw_subnet_size.append(10000)
        elif nw < 25: nw_subnet_size.append(13000)
        else: nw_subnet_size.append(16000)
subnet.remove_nodes_from(to_remove)

if lcc:
    subnet = max(nx.connected_component_subgraphs(subnet), key=len)
    addname = '_highnwLCC'

else:
    addname = '_highnw'
    degrees_zero = []
    for node in subnet.nodes:
        k = nx.degree(subnet, node)
        if k == 0:
            degrees_zero.append(node)
    subnet.remove_nodes_from(degrees_zero)
            
    
wij_list_sn = []
edge_colors_subnet = []
for u, v in subnet.edges:
    wij = subnet.get_edge_data(u, v)['weight']
    wij_list_sn.append(wij)
    if wij < 10: edge_colors_subnet.append('blue')
    elif wij < 20:  edge_colors_subnet.append('cyan')
    elif wij < 30:  edge_colors_subnet.append('greenyellow')
    elif wij < 40:  edge_colors_subnet.append('yellow')
    elif wij < 50:  edge_colors_subnet.append('orange')
    else: edge_colors_subnet.append('red')

plt.figure()
plt.hist(wij_list_sn, bins=range(1, max(wij_list) + 1, 3), normed=True)
#plt.scatter(xticks, hist)
plt.title('Pairwise weight distribution in the sub-network with $w/k\geq13$')
plt.xlabel('w(i,j)')
plt.ylabel('P(w(i, j))')
#plt.ylim(min(hist), max(hist))
plt.yscale('log')
plt.savefig('pairwise_distrib%s.pdf' %(addname))
plt.show() 

node_labels_subnet = {node: node_labels3[node] for node in subnet.nodes}

network_visualization.draw_network(subnet, node_labels_subnet, sizes=nw_subnet_size,
                                   edge_color_map=edge_colors_subnet,
                                   draw_labels=False,
                                   name='network_pictures/' + pdb_id +
                                   '_links_color' + addname, figsize=(50, 50),
                                   pos=pos_spring, facecolor='silver')


#subnetwork of wij > 40
subnet = net.copy()
lcc = True
to_remove = []
nw_subnet_size = []

for index, (u, v) in enumerate(subnet.edges):
    wij = subnet.get_edge_data(u, v)['weight']
    if wij < 30:
        to_remove.append((u, v))

subnet.remove_edges_from(to_remove)

if lcc:
    subnet = max(nx.connected_component_subgraphs(subnet), key=len)
    addname = '_highwijLCC'

else:
    addname = '_highwij'
    degrees_zero = []
    for node in subnet.nodes:
        k = nx.degree(subnet, node)
        if k == 0:
            degrees_zero.append(node)
    subnet.remove_nodes_from(degrees_zero)
            
for index, node in enumerate(subnet.nodes):
    k = nx.degree(subnet, node)
    w = nx.degree(subnet, node, weight='weight')
    nw = w / k
    if nw < 5: nw_subnet_size.append(1000)
    elif nw < 10: nw_subnet_size.append(4000)
    elif nw < 15: nw_subnet_size.append(7000)
    elif nw < 20: nw_subnet_size.append(10000)
    elif nw < 25: nw_subnet_size.append(13000)
    else: nw_subnet_size.append(16000)
        
wij_list_sn = []
edge_colors_subnet = []
for u, v in subnet.edges:
    wij = subnet.get_edge_data(u, v)['weight']
    wij_list_sn.append(wij)
    if wij < 40:  edge_colors_subnet.append('yellow')
    elif wij < 50:  edge_colors_subnet.append('orange')
    else: edge_colors_subnet.append('red')
    
plt.figure()
plt.hist(wij_list_sn, bins=range(1, max(wij_list) + 1), normed=True)
#plt.scatter(xticks, hist)
plt.title('Pairwise weight distribution in the sub-network with $w_{ij}\geq30$')
plt.xlabel('w(i,j)')
plt.ylabel('P(w(i, j))')
#plt.ylim(min(hist), max(hist))
plt.yscale('log')
plt.savefig('pairwise_distrib%s.pdf' %(addname))
plt.show() 

node_labels_subnet = {node: node_labels3[node] for node in subnet.nodes}

network_visualization.draw_network(subnet, node_labels_subnet,
                                   edge_color_map=edge_colors_subnet,
                                   draw_labels=True,
                                   name='network_pictures/' + pdb_id +
                                   '_links_color' + addname, figsize=(50, 50),
                                   pos=pos_spring, facecolor='silver')

#
#network_visualization.draw_network(subnet, node_labels3, sizes=nw_subnet_size,
#                                   edge_color_map=edge_colors_subnet,
#                                   draw_labels=False,
#                                   name='network_pictures/' + pdb_id +
#                                   '_links_color' + addname, figsize=(50, 50),
#                                   pos=pos_spring, facecolor='silver')