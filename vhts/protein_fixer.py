#!/usr/bin/env python

def read_pdb_protein(input_file, residue_index_initialize=True):

    fp = open(input_file)
    lines = fp.readlines()
    fp.close()

    protein_chain_residues = dict()

    pdb_info_lines = list()
#    ligand_info_lines = list()
    protein_dict = dict()
    ligand_dict = dict()
    conect_dict = dict()
    residue_num_old = ''
    residue_num2_old = ''
    chain_id_old = ''
    if residue_index_initialize:
        residue_num_new = 1
    for line in lines:
        if line[0:6] == 'HEADER':
            pdb_info_lines += [line]
            continue
        if line[0:6] == 'TITLE':
            pdb_info_lines += [line]
            continue
        if line[0:6] == 'DBREF ':
            pdb_info_lines += [line]
            continue
        if line[0:6] == 'SEQRES':
            pdb_info_lines += [line]
            continue

        if line[0:6] == 'ATOM  ':
            residue_name = line[17:20].strip()
            residue_num = int(line[22:26])
            residue_num2 = line[22:27]

            chain_id = line[21]
            altLoc = line[16]
            if altLoc != ' ' and altLoc != 'A':
                continue
            atom_number = int(line[6:11])
            if chain_id_old != chain_id:
                chain_id_old = chain_id
                residue_num_new = 1
                residue_num2_old = residue_num2
                residue_num_old = residue_num
            if residue_num2_old != residue_num2:
                residue_num_new += max(residue_num - residue_num_old, 1)
                residue_num2_old = residue_num2
                residue_num_old = residue_num

            line_out = '%s%4d %s' % (line[:22], residue_num_new, line[27:])
            protein_dict[atom_number] = line_out
            if chain_id not in protein_chain_residues:
                protein_chain_residues[chain_id] = dict()
            if residue_num2 not in protein_chain_residues[chain_id]:
                protein_chain_residues[chain_id][residue_num2] = residue_num_new

    for line in lines:
        if line[0:6] == 'HETATM':
            residue_name = line[17:20].strip()
            residue_num = int(line[22:26])
            residue_num2 = line[22:27]

            chain_id = line[21]
            atom_number = line[6:11]
            altLoc = line[16]
            if altLoc != ' ' and altLoc != 'A':
                continue
            if chain_id not in protein_chain_residues:
                continue
            if residue_num2 in protein_chain_residues[chain_id]:
                residue_num_new = protein_chain_residues[chain_id][residue_num2]
                line_out = '%s%4d %s' % (line[:22], residue_num_new, line[27:])
                protein_dict[atom_number] = 'ATOM  ' + line_out[6:]
            else:

            if chain_id_old != chain_id:
                chain_id_old = chain_id
                residue_num_new = 1
                residue_num2_old = residue_num2
                residue_num_old = residue_num
            if residue_num2_old != residue_num2:
                residue_num_new += max(residue_num - residue_num_old, 1)
                residue_num2_old = residue_num2
                residue_num_old = residue_num

                ligand_dict[atom_number] = line_out

    for line in lines:
        if line[0:6] == 'CONECT':
            conect_list = []
            for i in range(0, 8):
                ini = i * 5 + 6
                fin = (i + 1) * 5 + 6
                atom_number = line[ini:fin].strip()
                if len(atom_number) > 0:
                    conect_list += [int(atom_number)]
            conect_idx = conect_list[0]
            if conect_idx not in conect_dict:
                conect_dict[conect_idx] = conect_list[1:]
            else:
                conect_dict[conect_idx] = conect_dict[conect_idx] + \
                    conect_list[1:]

    return (protein_chain_residues, pdb_info_lines, protein_dict, conect_dict)

def separate_pdb(input_file, protein_file):
    """
        This function is for residue_re-numbering
        example A:4, A:4a, A:5, ...
            -> A:1, A:2, A:3, ...
        pdbfixer does not identify 4, 4a...
        input:
            input pdb file
            output pdb file
        return:
            None
    """
    result = read_pdb_protein(input_file, residue_index_initialize=True)
    (protein_chain_residues, pdb_info_lines, protein_dict, conect_dict) = result

    protein_atom_numbers = sorted(protein_dict.keys())
    fp_p = open(protein_file, 'w')
    for line in pdb_info_lines:
        line_out = line
        fp_p.write(line_out)

    for atom_number in protein_atom_numbers:
        line_out = protein_dict[atom_number]
        fp_p.write(line_out)

    for atom_number in protein_atom_numbers:
        if atom_number not in conect_dict:
            continue
        ans = conect_dict[atom_number]
        ans2 = list()
        for an in ans:
            if an in protein_atom_numbers:
                ans2 += [an]
        num_conect = len(ans2)
        line_out = ''
        for i_con in range(num_conect):
            if i_con % 4 == 0:
                line_out += 'CONECT%5d' % (atom_number)
            line_out += '%5d' % (ans2[i_con])
            if i_con % 4 == 3:
                line_out += '\n'
        if line_out[-1] != '\n':
            line_out += '\n'
#        if len(line_out) < 14:
#            continue
        fp_p.write(line_out)
    line_out = 'END'
    fp_p.write(line_out)

    fp_p.close()

def fix_pdb(protein_file, protein_fix_file, pH=7.4):

    fixer = pdbfixer.PDBFixer(filename=protein_file)
#    fixer.applyMutations(["ALA-57-GLY"], "A")
#    fixer.missingResidues = {}
    fixer.findMissingResidues()
    missing_residues = fixer.missingResidues

    fixer.findNonstandardResidues()
    nonstandard_residues = fixer.nonstandardResidues
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    missing_atoms = fixer.missingAtoms
    missing_terminals = fixer.missingTerminals
    fixer.addMissingAtoms()
    try:
        fixer.addMissingHydrogens(pH=pH)
    except:
        print('Error: fixer.addMissingHydrogens(pH=pH)', protein_file)
    fp_out = open(protein_fix_file, 'w')
    app.PDBFile.writeFile(fixer.topology, fixer.positions, fp_out, True)
    fp_out.close()

    return (missing_residues, nonstandard_residues,
            missing_atoms, missing_terminals)

def pybel_rewrite_pdb(input_file, output_file, mol_format='pdb'):
    ms = pybel.readfile(mol_format, input_file)
    m = list(ms)[0]
    m.write('pdb', output_file, overwrite=True)


def fix_charge_pdb(protein_conect_file, protein_charge_file):
    result = read_pdb_protein(
        protein_conect_file, residue_index_initialize=False)
    (protein_chain_residues, pdb_info_lines, protein_dict, conect_dict) = result

    protein_atom_numbers = sorted(protein_dict.keys())
    fp_p = open(protein_charge_file, 'w')
    for line in pdb_info_lines:
        line_out = line
        fp_p.write(line_out)

    for atom_number in protein_atom_numbers:
        line_out = protein_dict[atom_number]
        if atom_number not in conect_dict:
#            line_out2 = line_out
#            fp_p.write(line_out2)
            continue
        ans = conect_dict[atom_number]
        ans2 = list()
        for an in ans:
            if an in protein_atom_numbers:
                ans2 += [an]
        num_conect = len(ans2)
        atom_type = line_out[77]
        charge = 0
        if atom_type == 'N':
            charge = num_conect - 3
        if atom_type == 'O':
            charge = num_conect - 2

        if charge > 0:
            line_out2 = '%s%1d+\n' % (line_out[0:78], charge)
        elif charge < 0:
            line_out2 = '%s%1d-\n' % (line_out[0:78], -charge)
        else:
            line_out2 = line_out
        fp_p.write(line_out2)
    for atom_number in protein_atom_numbers:
        if atom_number not in conect_dict:
            continue
        ans = conect_dict[atom_number]
        ans2 = list()
        for an in ans:
            if an in protein_atom_numbers:
                ans2 += [an]
        num_conect = len(ans2)
        line_out = ''
        for i_con in range(num_conect):
            if i_con % 4 == 0:
                line_out += 'CONECT%5d' % (atom_number)
            line_out += '%5d' % (ans2[i_con])
            if i_con % 4 == 3:
                line_out += '\n'
        if line_out[-1] != '\n':
            line_out += '\n'
#        if len(line_out) < 14:
#            continue
        fp_p.write(line_out)
    line_out = 'END'
    fp_p.write(line_out)
    fp_p.close()


def main():

    import argparse
    title_line = 'Fixer for protein pdb which is converted from pdbqt'
    parser = argparse.ArgumentParser(description=title_line)
    parser.add_argument('-i', '--input_file', required=True,
                        help='input protein pdb file')
    parser.add_argument('-o', '--output_file', required=True,
                        help='output protein pdb file')
    args = parser.parse_args()
    protein_input_file = args.input_file
    protein_output_file = args.output_file

    separate_pdb(protein_input_file, protein_output_file)
    result = fix_pdb(protein_output_file, protein_output_file, pH=7.4)
    (missing_residues, nonstandard_residues, missing_atoms, missing_terminals) = result
    pybel_rewrite_pdb(protein_output_file, protein_output_file)
    fix_charge_pdb(protein_output_file, protein_output_file)


if __name__ == "__main__":
    main()
