import biographs as bg
from biopandas.pdb import PandasPdb
import numpy as np

def assign_secondary_structure(pdb):
    
    ppdb = PandasPdb().read_pdb(pdb)
    
    secondary_structure = {}
    
    helices_from_pdb = ppdb.df['OTHERS'][ppdb.df['OTHERS']['record_name'] == 'HELIX']['entry']
    for helix in helices_from_pdb:
        identifier_h = helix[5:8].strip()
        initial_chain_h = helix[13].strip()
        initial_pos_h = helix[16:19].strip()
        final_pos_h = helix[28:31].strip()
        for i in range(int(initial_pos_h), int(final_pos_h) + 1):
            secondary_structure[initial_chain_h + str(i)] = 'helix' + identifier_h + '-' + initial_chain_h
            
    sheets_from_pdb = ppdb.df['OTHERS'][ppdb.df['OTHERS']['record_name'] == 'SHEET']['entry']
    for sheet in sheets_from_pdb:
        identifier_s = sheet[6:8].strip()
        initial_chain_s = sheet[15].strip()
        initial_pos_s = sheet[17:20].strip()
        final_pos_s = sheet[28:31].strip()
        for i in range(int(initial_pos_s), int(final_pos_s) + 1):
            secondary_structure[initial_chain_s + str(i)] = 'sheet' + identifier_s + '-' + initial_chain_s

    mol = bg.Pmolecule(pdb)
    net = mol.network()
    
    residues_type = {}
    for residue in mol.model.get_residues():
        res_type = residue.resname
        res_pos = residue.parent.id + str(residue.id[1])
        residues_type[res_pos] = res_type
    
    residues = list(net.nodes) #assume they are ordered
    last_structure = None
    last_chain = None
    i = 0
    for residue in residues:
        chain = residue[0]
        try:
            structure = secondary_structure[residue]
            if structure != last_structure:
                i += 1
        except KeyError:
            if chain != last_chain:
                i += 1
            structure = 'loop' + str(i)
            secondary_structure[residue] = structure
        last_structure = structure
        last_chain = chain


    return secondary_structure

def get_neighbor_structure_relation(secondary_structure, u, v):
    
    chain_u = u[0]
    chain_v = v[0]
    pos_u = u[1::]
    pos_v = v[1::]
    struct_u = secondary_structure[u]
    struct_v = secondary_structure[v]
    
    if chain_u == chain_v:
        if struct_u == struct_v:
            dist = np.abs(int(pos_u) - int(pos_v))
            if dist < 5:
                relation = '2D' + str(dist)
            else:
                relation = '3D'
        else:
            relation = '3D'
    else:
        relation = '4D'
    
    return relation