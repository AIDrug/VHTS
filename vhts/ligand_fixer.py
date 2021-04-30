#!/usr/bin/env python


def read_pdb_ligand(input_file):

    fp = open(input_file)
    lines = fp.readlines()
    fp.close()

    model_dict = dict()

    pdb_info_lines = list()
    ligand_dict = dict()
    conect_dict = dict()
    for line in lines:
        if line[0:6] == 'MODEL ':
            model = int(line[6:].strip())
            pdb_info_lines = list()
            ligand_dict = dict()
            conect_dict = dict()

        if line[0:6] == 'HEADER':
            pdb_info_lines += [line]
            continue
        if line[0:6] == 'TITLE':
            pdb_info_lines += [line]
            continue
        if line[0:6] == 'REMARK':
            pdb_info_lines += [line]
            continue
        if line[0:6] == 'COMPND':
            pdb_info_lines += [line]
            continue
        if line[0:6] == 'AUTHOR':
            pdb_info_lines += [line]
            continue

        if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
            altLoc = line[16]
            if altLoc != ' ' and altLoc != 'A':
                continue
            atom_number = int(line[6:11])
            line_out = line
            ligand_dict[atom_number] = line_out

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
        if line[0:6] == 'ENDMDL':
            model_dict[model] = (pdb_info_lines, ligand_dict, conect_dict)
    if len(model_dict.keys()) == 0:
        model_dict[1] = (pdb_info_lines, ligand_dict, conect_dict)

    return (model_dict)


def fix_one_model(ligand_model):
    (pdb_info_lines, ligand_dict, conect_dict) = ligand_model
    line_out_total = ''
    ligand_atom_numbers = sorted(ligand_dict.keys())
    for line in pdb_info_lines:
        line_out = line
        line_out_total += line_out

    for atom_number in ligand_atom_numbers:
        line = ligand_dict[atom_number]
        if atom_number not in conect_dict:
            continue
        atom_type = line[77]
        if atom_type == 'N':
            ans = conect_dict[atom_number]
            con_dict = dict()
            for an in ans:
                if an in ligand_atom_numbers:
                    if an not in con_dict:
                        con_dict[an] = 0
                    con_dict[an] += 1
            con_dict2 = dict()
            ans2 = list()
            keys = sorted(con_dict.keys())
            check_o2 = False
            check_mod_o2 = False
            o2_atom_num = 0
            for key in keys:
                bond_type = con_dict[key]
                line2 = ligand_dict[key]
                atom_type2 = line2[77]
                if atom_type2 == 'O' and bond_type == 2 and not check_o2:
                    check_o2 = True
                elif atom_type2 == 'O' and bond_type == 2 and check_o2:
                    bond_type = 1
                    o2_atom_num = key
                    check_mod_o2 = True
                con_dict2[key] = bond_type
                for i in range(bond_type):
                    ans2 += [key]
            conect_dict[atom_number] = ans2
            if check_mod_o2:
                ans_o = conect_dict[o2_atom_num]
                con_dict = dict()
                for an in ans_o:
                    if an in ligand_atom_numbers:
                        if an not in con_dict:
                            con_dict[an] = 0
                        con_dict[an] += 1
                ans_o2 = list()
                keys = sorted(con_dict.keys())
                for key in keys:
                    bond_type = con_dict[key]
                    if key == atom_number and bond_type == 2:
                        bond_type = 1
                    for i in range(bond_type):
                        ans_o2 += [key]
                conect_dict[o2_atom_num] = ans_o2

    for atom_number in ligand_atom_numbers:
        line_out = ligand_dict[atom_number]
        if atom_number not in conect_dict:
            continue
        ans = conect_dict[atom_number]
        ans2 = list()
        for an in ans:
            if an in ligand_atom_numbers:
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
        line_out_total += line_out2

    for atom_number in ligand_atom_numbers:
        if atom_number not in conect_dict:
            continue
        ans = conect_dict[atom_number]
        ans2 = list()
        for an in ans:
            if an in ligand_atom_numbers:
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
#        if len(line_out) < 1
#            continue
        line_out_total += line_out
    line_out = 'END\n'
    line_out_total += line_out
    return line_out_total


def fix_ligand_pdb(ligand_input_file, ligand_fix_file):
    model_dict = read_pdb_ligand(ligand_input_file)
    fp_p = open(ligand_fix_file, 'w')
    keys = sorted(model_dict.keys())
    num_model = len(keys)
    for model in keys:
        ligand_model = model_dict[model]
        if num_model > 1:
            line_out = 'MODEL %8d\n' % model
            fp_p.write(line_out)
        line_out = fix_one_model(ligand_model)
        fp_p.write(line_out)
        if num_model > 1:
            line_out = 'ENDMDL\n'
            fp_p.write(line_out)

    fp_p.close()


def pdbqt_to_pdb(input_pdbqt_file, output_pdb_file):
    run_line = 'obabel %s -O %s' % (input_pdbqt_file, output_pdb_file)
    result = subprocess.check_output(run_line.split(),
                                     stderr=subprocess.STDOUT,
                                     universal_newlines=True)
    fix_ligand_pdb(output_pdb_file, output_pdb_file)
    run_line = 'obabel %s -h -O %s' % (output_pdb_file, output_pdb_file)
    result = subprocess.check_output(run_line.split(),
                                     stderr=subprocess.STDOUT,
                                     universal_newlines=True)
    return

def main():

    import argparse
    title_line = 'Fixer for ligand pdb which is converted from pdbqt'
    parser = argparse.ArgumentParser(description=title_line)
    parser.add_argument('-i', '--input_file', required=True,
                        help='input ligand pdb file')
    parser.add_argument('-o', '--output_file', required=True,
                        help='output ligand pdb file')
    args = parser.parse_args()
    ligand_input_file = args.input_file
    ligand_output_file = args.output_file

    fix_ligand_pdb(ligand_input_file, ligand_output_file)


if __name__ == "__main__":
    main()
