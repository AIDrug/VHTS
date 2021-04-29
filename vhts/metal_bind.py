#!/usr/bin/env python
import argparse
import numpy as np
from openbabel import pybel
import pandas as pd


def get_bond_info(m):
    bond_dict = dict()
    mol = m.OBMol

    for i in range(mol.NumBonds()):
        bb = mol.GetBondById(i)
        if bb is None:
            continue
        begin = bb.GetBeginAtomIdx()
        end = bb.GetEndAtomIdx()
        if begin not in bond_dict:
            bond_dict[begin] = [end]
        else:
            bond_dict[begin] += [end]
        if end not in bond_dict:
            bond_dict[end] = [begin]
        else:
            bond_dict[end] += [begin]
    return bond_dict


def find_neighbor_metal(ligand_file, metal_coor, dist_cutoff=2.5, skip_neighbor_hydrogen=True):

    ms = pybel.readfile('pdb', ligand_file)
    ms= list(ms)
    bond_dict = get_bond_info(ms[0])

    num_conformers = len(ms)
    if num_conformers>1:
        num_conformers = num_conformers - 1

    num_bind_atom_list = list()
    num_bind_max= 0
    model_num = 0

    for i in range(num_conformers):
        m = ms[i]
        atoms = m.atoms
        bind_list = list()
        for atom in atoms:
            atom_idx = atom.idx
            atomic_num = atom.atomicnum
            neighbor_hydrogen = False
            if atomic_num==7 or atomic_num==8:
#                if skip_neighbor_hydrogen:
                if skip_neighbor_hydrogen and atomic_num==7:
                    neighbor_atoms_idx = bond_dict[atom_idx]
                    for neighbor_atom_idx in neighbor_atoms_idx:
                        neighbor_atom = atoms[neighbor_atom_idx - 1]
                        if neighbor_atom.atomicnum == 1:
                            neighbor_hydrogen = True
                    if neighbor_hydrogen :
                        continue
                coor = np.array(atom.coords)
                dist = np.linalg.norm(metal_coor-coor)
                if dist < dist_cutoff:
                    bind_list += [[atom_idx, atomic_num, dist]]

        num_bind_atom = 0
        batom_idx = list()
        for bind in bind_list:
            atom_idx, atomic_num, dist = bind
            atom = atoms[atom_idx-1]
            neighbor_atoms_idx = bond_dict[atom_idx]
            neighbor_check = False
            for neighbor_atom_idx in neighbor_atoms_idx:
                if neighbor_atom_idx in batom_idx:
                    neighbor_check = True
            if not neighbor_check :
                num_bind_atom += 1
            batom_idx += [atom_idx]
        num_bind_atom_list += [num_bind_atom]
    return num_bind_atom_list


def main():

#    ligand_file = 'pdb/2G1MA_4HG.pdb'

    metal_coor = np.array([39.864, 22.370, 13.612])
#    ligand_list = ['2G1MA_4HG', '4JZRA_4JR', '4KBZA_1QA',
#                    '5V18A_8UY', '6ST3A_LUW', '6YVTA_PW2']

    df = pd.read_csv('docking.csv')
    num_data = df.shape[0]
    size = (400,400)
    num_d = 100
    df2 = df.sort_values('Docking1')[0:num_d]
    for i in range(num_d):
        mdata = df2.iloc[i]
        ligand = mdata['MOL_ID']
        docking1 = mdata['Docking1']
        docking = mdata['Docking'].strip('[').strip(']').split(',')
        dir1 = ligand[0:4]
        ligand_file = 'docking/%s/dock_%s.pdb' % (dir1, ligand)

        num_bind_atom_list = find_neighbor_metal(ligand_file, metal_coor)
        line_out = '%s %s %7.3f %s' % (ligand, num_bind_atom_list, docking1, docking)
        print(line_out)

if __name__ == '__main__':
    main()

