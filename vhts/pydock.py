#!/usr/bin/env python
import sys
import os
import numpy as np
from multiprocessing import Manager
from multiprocessing import Process
from multiprocessing import Queue
import subprocess
from openbabel import pybel
import argparse
import pandas as pd
from vhts import metal_bind
from vhts import ligand_tools


class DockingVina(object):
    """
    python module for Vina
    """

    def __init__(self, docking_params):
        """
            Construction Docking object
            put parameters with dictionary
            docking_params['vina_program']
            docking_params['dock_config_file']
            docking_params['output_save']
            docking_params['gen3d_dir']
            docking_params['dock_dir']
            docking_params['num_sub_proc']
            docking_params['timeout_gen3d']
            docking_params['timeout_dock']
            docking_params['tlen']
            docking_params['pout']

        """
        self.vina_program = docking_params['vina_program']

        self.num_sub_proc = docking_params['num_sub_proc']
        self.timeout_gen3d = docking_params['timeout_gen3d']
        self.timeout_dock = docking_params['timeout_dock']
        self.neutralize = docking_params['neutralize']
        self.pH = docking_params['pH']
        self.dock_config_file = docking_params['dock_config_file']

        self.output_save = docking_params['output_save']
        if self.output_save:
            self.tlen = docking_params['tlen']

        self.pout = docking_params['pout']

        self.gen3d_dir = docking_params['gen3d_dir']
        if not os.path.exists(self.gen3d_dir):
            try:
                os.makedirs(self.gen3d_dir)
            except FileExistsError as e:
                print(e, flush=True)

        self.dock_dir = docking_params['dock_dir']
        if not os.path.exists(self.dock_dir):
            try:
                os.makedirs(self.dock_dir)
            except FileExistsError as e:
                print(e, flush=True)

        self.use_my_module = docking_params['use_my_module']
        if self.use_my_module:
            my_module_path = docking_params['my_module_path']
            my_module_dir = os.path.dirname(my_module_path)
            sys.path.append(my_module_dir)
            import my_module
            self.my_class = my_module.my_class
            self.my_class.__init__(self, docking_params)

        self.check_metal_bind = docking_params['check_metal_bind']
        if self.check_metal_bind:
            self.metal_coor = docking_params['metal_coor']
            self.metal_bind_cutoff = docking_params['metal_bind_cutoff']
        self.rescoring = docking_params['rescoring']
        self.rescoring_program = docking_params['rescoring_program']
        self.rescoring_config_file = docking_params['rescoring_config_file']

    def docking(self, ligand_file, docking_pdbqt_file, docking_log_file):
        """
            run_docking program using subprocess
            input :
                ligand_file
                docking_pdbqt_file
            output :
                affinity list for a input molecule
        """

        run_line = '%s' % self.vina_program
        run_line += ' --config %s' % self.dock_config_file
        run_line += ' --ligand %s' % ligand_file
        run_line += ' --out %s' % docking_pdbqt_file
        if self.output_save:
            run_line += ' --log %s' % (docking_log_file)
        e = None
        try:
            result = subprocess.check_output(run_line.split(),
                                             stderr=subprocess.STDOUT,
                                             timeout=self.timeout_dock,
                                             universal_newlines=True)
        except Exception as e:
            return [99.999], e

        result_lines = result.split('\n')

        check_result = False
        affinity_list = list()
        for result_line in result_lines:
            if result_line.startswith('-----+'):
                check_result = True
                continue
            if not check_result:
                continue
            if result_line.startswith('Writing output'):
                break
            if result_line.startswith('Refine time'):
                break
            lis = result_line.strip().split()
            if not lis[0].isdigit():
                break
#            mode = int(lis[0])
            affinity = float(lis[1])
            affinity_list += [affinity]
        if len(affinity_list) == 0:
            e = 'WARNING: Could not find any conformations.'
            return [99.999], e
        return affinity_list, e

    def docking_score_only(self, ligand_file):
        """
            run docking program with score_only using subprocess
            input :
                ligand_file
            output :
                affinity list for a input molecule
        """

        run_line = '%s' % self.rescoring_program
        run_line += ' --config %s' % self.rescoring_config_file
        run_line += ' --ligand %s' % ligand_file
        run_line += ' --score_only'

        e = None
        try:
            result = subprocess.check_output(run_line.split(),
                                             stderr=subprocess.STDOUT,
                                             timeout=self.timeout_dock,
                                             universal_newlines=True)
        except Exception as e:
            return [99.999], e

        result_lines = result.split('\n')

#        weight_list = list()
#        check_weight = False
        affinity_list = list()
        for result_line in result_lines:
            #            if result_line.startswith('Weights'):
            #                check_weight = True
            #                continue
            #            if check_weight:
            #                lis = result_line.strip().split()
            #                if len(lis) <2:
            #                    check_weight = False
            #                    continue
            #                weight_list += [[float(lis[0]), lis[1]]]
            #                continue
            if result_line.startswith('Affinity:'):
                lis = result_line.strip().split()
                affinity = float(lis[1])
                affinity_list += [affinity]
        if len(affinity_list) == 0:
            return [99.999], e
        return affinity_list, e

    def creator(self, q, data, num_sub_proc):
        """
            put data to queue
            input: queue
                data = [(idx1, molid1, smi1), (idx2, molid2, smi2), ...]
                num_sub_proc (for end signal)
        """
        for d in data:
            idx = d[0]
            q.put((idx, d[1]))

        for i in range(0, num_sub_proc):
            q.put('DONE')

    def simulation_process(self, idx, mol_id, smi, pid):

        result_dict = dict()
        if self.neutralize or (self.pH is not None):
            smi_p = ligand_tools.ligand_preparation(smi, self.neutralize, self.pH)
        else:
            smi_p = smi
        if not self.output_save:
            ligand_pdb_file = '%s/ligand_%d.pdb' % (self.gen3d_dir, pid)
            ligand_pdbqt_file = '%s/ligand_%s.pdbqt' % (self.gen3d_dir, pid)
            docking_pdbqt_file = '%s/dock_%d.pdbqt' % (
                self.dock_dir, pid)
            docking_log_file = '%s/dock_%d.log' % (self.dock_dir, pid)
            docking_pdb_file = '%s/dock_%s.pdb' % (self.dock_dir, pid)
            out_dock_dir1 = None
        else:
            mol_id2 = mol_id[0:self.tlen]
            out_gen3d_dir1 = self.gen3d_dir + "/" + mol_id2
            if not os.path.exists(out_gen3d_dir1):
                try:
                    os.makedirs(out_gen3d_dir1)
                except FileExistsError as e:
                    print(e, flush=True)
            ligand_pdb_file = '%s/ligand_%s.pdb' % (out_gen3d_dir1, mol_id)
            ligand_pdbqt_file = '%s/ligand_%s.pdbqt' % (out_gen3d_dir1, mol_id)

            out_dock_dir1 = self.dock_dir + "/" + mol_id2
            if not os.path.exists(out_dock_dir1):
                try:
                    os.makedirs(out_dock_dir1)
                except FileExistsError as e:
                    print(e, flush=True)

            docking_pdbqt_file = '%s/dock_%s.pdbqt' % (out_dock_dir1, mol_id)
            docking_pdb_file = '%s/dock_%s.pdb' % (out_dock_dir1, mol_id)
            docking_log_file = '%s/dock_%s.log' % (out_dock_dir1, mol_id)

        e = ligand_tools.gen_3d(smi_p, ligand_pdb_file)
        if e is not None:
            e2 = ligand_tools.gen_3d(smi_p, ligand_pdb_file)
            if e2 is not None:
                print(e2, 'gen_3d', idx, mol_id, smi_p, flush=True)
                docking_score = np.array([99.999], dtype=np.float32)
                result_dict['docking'] = docking_score
                return result_dict
        e = ligand_tools.pdb_to_pdbqt(ligand_pdb_file, ligand_pdbqt_file)
        if e is not None:
            print(e, 'pdb_to_pdbqt', idx, mol_id, smi_p, flush=True)
            docking_score = np.array([99.999], dtype=np.float32)
            result_dict['docking'] = docking_score
            return result_dict

        docking_score, e = self.docking(ligand_pdbqt_file, docking_pdbqt_file,
                                        docking_log_file)
        docking_score = np.array(docking_score, dtype=np.float32)
        if e is not None:
            docking_score = [99.999]
            result_dict['docking'] = docking_score
            print(e, 'docking', idx, mol_id, smi_p, flush=True)
            return result_dict
        result_dict['docking'] = docking_score
        if self.output_save or self.check_metal_bind or self.rescoring or self.use_my_module:
            ligand_tools.pdbqt_to_pdb_ref(
                docking_pdbqt_file, docking_pdb_file, ligand_pdb_file)

        if self.rescoring:
            docking_rescore, e = self.docking_score_only(docking_pdb_file)
            docking_rescore = np.array(docking_rescore, dtype=np.float32)
            if e is not None:
                docking_rescore = np.array([99.999], dtype=np.float32)
                print(e, 're-scoring', idx, mol_id, smi_p, flush=True)
            result_dict['docking_re'] = docking_rescore

        if self.use_my_module:
            self.my_class.simulation_process(self, idx, mol_id, smi, smi_p, pid,
                                             out_dock_dir1, docking_pdb_file,
                                             result_dict)

        if self.check_metal_bind:
            try:
                result = metal_bind.find_neighbor_metal(
                    docking_pdb_file, self.metal_coor,
                    dist_cutoff=self.metal_bind_cutoff,
                    skip_neighbor_hydrogen=True)
                num_metal_bind_atom_list = np.array(result, dtype=np.float32)
            except Exception as e:
                num_metal_bind_atom_list = np.array([0], dtype=np.float32)
                print(e, 'check metal bind', idx, mol_id, smi_p, flush=True)
            result_dict['num_metal_bind_atom'] = num_metal_bind_atom_list

        return result_dict

    def worker(self, q, return_dict):
        """
            generate subprocess for docking
            input
                q (queue)
                return_dict
        """
        pid = os.getpid()
        while True:
            qqq = q.get()
            if qqq == 'DONE':
                # print('proc =', os.getpid())
                break

            (idx, d) = qqq
            mol_id = d[0]
            smi = d[1]
            # print screening processing in every pout step
            if self.pout != 0:
                if idx % self.pout == self.pout-1:
                    print("processing: ", idx+1, flush=True)
            result_dict = self.simulation_process(idx, mol_id, smi, pid)
            return_dict[idx] = result_dict

    def predict(self, smiles_list):
        """
            input SMILES list
            output result_dict
            result_dict include affinity list (and other scores)
            corresponding to the SMILES list
            if docking is fail, docking score is [99.999]
        """
        data = list(enumerate(smiles_list))
        num_data = len(data)
        num_sub_proc = min(self.num_sub_proc, num_data)

        q1 = Queue()
        manager = Manager()
        return_dict = manager.dict()
        proc_master = Process(target=self.creator,
                              args=(q1, data, num_sub_proc))
        proc_master.start()

        # create slave process
        procs = []
        for sub_id in range(0, num_sub_proc):
            proc = Process(target=self.worker, args=(q1, return_dict))
            procs.append(proc)
            proc.start()

        q1.close()
        q1.join_thread()
        proc_master.join()
        for proc in procs:
            proc.join()
        keys = sorted(return_dict.keys())

        result_dict = dict()
        docking_score_list = list()
        if self.check_metal_bind:
            num_metal_bind_atom_list = list()
        if self.rescoring:
            docking_re_list = list()

        for key in range(num_data):
            if key in keys:
                result_dict0 = return_dict[key]
                if 'docking' in result_dict0:
                    docking_score = result_dict0['docking']
                else:
                    docking_score = np.array([99.999], dtype=np.float32)

                if self.check_metal_bind:
                    if 'num_metal_bind_atom' in result_dict0:
                        num_metal_bind_atom = result_dict0['num_metal_bind_atom']
                    else:
                        num_metal_bind_atom = np.array([0], dtype=np.float32)
                if self.rescoring:
                    if 'docking_re' in result_dict0:
                        docking_re = result_dict0['docking_re']
                    else:
                        docking_re = np.array([99.999], dtype=np.float32)

            else:
                docking_score = np.array([99.999], dtype=np.float32)
                if self.check_metal_bind:
                    num_metal_bind_atom = np.array([0], dtype=np.float32)
                if self.rescoring:
                    docking_re = np.array([99.999], dtype=np.float32)

            docking_score_list += [docking_score]
            if self.check_metal_bind:
                num_metal_bind_atom_list += [num_metal_bind_atom]
            if self.rescoring:
                docking_re_list += [docking_re]

        result_dict['docking'] = docking_score_list
        if self.check_metal_bind:
            result_dict['num_metal_bind_atom'] = num_metal_bind_atom_list
        if self.rescoring:
            result_dict['docking_re'] = docking_re_list

        if self.use_my_module:
            self.my_class.predict(self, smiles_list, result_dict, return_dict)

        return result_dict


def cal_box_size(ligand_file_list, margin=4.0, use_hydrogen=False):
    """
        cal box size from ligands
        input:
            ligand file list
            margin: addtional box-size to ligand size, default 3.0
            use_hydrogen: include hydrogen atom position, defalut Flase
        output:
            box_center : tuple (x, y, z)
            box_size : tuple (wx, wy, wz)
    """
    cmins = list()
    cmaxs = list()
    for ligand_file in ligand_file_list:
        file_format = ligand_file.split(".")[-1]
        ms = list(pybel.readfile(file_format, ligand_file))
        m = ms[0]
        if not use_hydrogen:
            m.removeh()
        atoms = m.atoms
        coor_list = list()
        for atom in atoms:
            coor_list.append(atom.coords)
        coor = np.array(coor_list)

        cmin0 = coor.min(axis=0)
        cmax0 = coor.max(axis=0)
        cmins.append(cmin0)
        cmaxs.append(cmax0)
    cmins = np.array(cmins)
    cmaxs = np.array(cmaxs)
    cmin = cmins.min(axis=0)
    cmax = cmaxs.max(axis=0)

    box_center = tuple((cmax+cmin)/2.0)
    box_size = tuple((cmax-cmin) + margin*2)

    return box_center, box_size


class LoadFromConfig(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        with values as f:
            parser.parse_args(f.read().split(), namespace)


class ExtendAction(argparse.Action):

    def __call__(self, parser, namespace, values, option_string=None):
        items = getattr(namespace, self.dest) or []
        items.extend(values)
        setattr(namespace, self.dest, items)


def parser_arg(parser):
    # vina parameter

    parser.register('action', 'extend', ExtendAction)
    parser.add_argument('--arg_file', type=open, required=False, default=None,
                        action=LoadFromConfig, help='argment file')
    parser.add_argument('--dock_config', type=str, required=False,
                        default=None, help='docking config file ')
    parser.add_argument('-v', '--vina_program', type=str, required=False,
                        default='qvina02',
                        help='select vina, qvina02, or smina')
    parser.add_argument('--my_module', type=str, required=False,
                        default=None,
                        help='set user python module path (for pifinder)')
    parser.add_argument('--neutralize', action='store_true',
                        required=False, help='neutralize smiles ')
    parser.add_argument('--pH', type=float, default=None,
                        required=False, help='protonate state for pH 7.4 ')
    parser.add_argument('--output_save', action='store_true', required=False,
                        help='default output pdbqt is temp file ')
    parser.add_argument('--gen3d_dir', type=str, required=False, default='tmp',
                        help='3d initial conformation directory')
    parser.add_argument('--dock_dir', type=str, required=False,
                        default='tmp', help='binding conformation directory')
    parser.add_argument('--num_sub_proc', type=int, required=False,
                        default=10, help=' --num_sub_proc 10')
    parser.add_argument('--timeout_gen3d', type=int, required=False,
                        default=1, help=' --timeout_gen3d 1')
    parser.add_argument('--timeout_dock', type=int, required=False,
                        default=120, help=' --timeout_dock 120')
    parser.add_argument('--tlen', type=int, default='7', required=False,
                        help='lenth of sub directory name, default: 7')
    parser.add_argument('--pout', type=int, default='0', required=False,
                        help='print processing out: 0 or number, default: 0')
    parser.add_argument('--metal_coor', type=str, default=None, required=False,
                        help='position of metal_ion,' +
                        ' example: --metal_coor="1.0,-1.0,0.0" default: None')
    parser.add_argument('--metal_cutoff', type=float, required=False,
                        default=3.0,
                        help='metal ion - HDA cutoff distance,')
    parser.add_argument('--rescoring_program', type=str, required=False,
                        default='smina', help='smina path')
    parser.add_argument('--rescoring_config', type=str, required=False,
                        default=None, help='docking config file for rescoring')

    return


def arg_to_params(parser):

    use_my_module = False
    for i, m in enumerate(sys.argv):
        if m == '--my_module':
            my_module_path = sys.argv[i+1]
            use_my_module = True
            my_module_dir = os.path.dirname(my_module_path)
            sys.path.append(my_module_dir)
            import my_module
            my_module.parser_arg(parser)

    args = parser.parse_args()

    vina_program = args.vina_program
    num_sub_proc = args.num_sub_proc
    timeout_gen3d = args.timeout_gen3d
    timeout_dock = args.timeout_dock
    output_save = args.output_save
    gen3d_dir = args.gen3d_dir
    dock_dir = args.dock_dir
    dock_config_file = args.dock_config

    tlen = args.tlen
    pout = args.pout
    neutralize = args.neutralize
    pH = args.pH

    check_metal_bind = False
    if args.metal_coor is not None:
        check_metal_bind = True
        metal_coor = np.array(args.metal_coor.strip(
            '"').split(','), dtype=np.float32)
        if metal_coor.shape[0] != 3:
            print('metal coordinate is strange', args.metal_coor)
            sys.exit()
        metal_bind_cutoff = args.metal_cutoff

    rescoring = False
    rescoring_config_file = args.rescoring_config
    rescoring_program = args.rescoring_program
    if rescoring_config_file is not None:
        rescoring = True

    docking_params = dict()
    docking_params['vina_program'] = vina_program
    docking_params['gen3d_dir'] = gen3d_dir
    docking_params['dock_dir'] = dock_dir
    docking_params['num_sub_proc'] = num_sub_proc
    docking_params['timeout_gen3d'] = timeout_gen3d
    docking_params['timeout_dock'] = timeout_dock
    docking_params['output_save'] = output_save
    docking_params['tlen'] = tlen
    docking_params['pout'] = pout
    docking_params['neutralize'] = neutralize
    docking_params['pH'] = pH
    docking_params['dock_config_file'] = dock_config_file
    docking_params['check_metal_bind'] = check_metal_bind
    if args.metal_coor is not None:
        docking_params['metal_coor'] = metal_coor
        docking_params['metal_bind_cutoff'] = metal_bind_cutoff
    docking_params['rescoring'] = rescoring
    docking_params['rescoring_program'] = rescoring_program
    docking_params['rescoring_config_file'] = rescoring_config_file

    my_module_path = args.my_module
    docking_params['use_my_module'] = use_my_module
    docking_params['my_module_path'] = my_module_path

    if use_my_module:
        docking_params = my_module.arg_to_params(parser, docking_params)

    return args, docking_params


def main():

    parser = argparse.ArgumentParser(description='docking with multi process')
    parser.add_argument('-l', '--ligand_list_file', type=str, required=False,
                        default='smiles.txt',
                        help=' --ligand_list_file smiles.txt')
    parser.add_argument('-o', '--output_file', type=str, required=False,
                        default='docking.txt',
                        help=' --output_file docking.txt')

    parser_arg(parser)

    args, docking_params = arg_to_params(parser)
    ligand_list_file = args.ligand_list_file
    output_file = args.output_file

    if len(sys.argv) < 2:
        parser.print_usage()
        sys.exit()
    if args.dock_config is None:
        parser.print_usage()
        print('dock_config is missing')
        sys.exit()

    ligand_file_format = ligand_list_file.strip().split('.')[-1]
    if ligand_file_format == 'txt':
        field_separator = '\s+'
    elif ligand_file_format == 'csv':
        field_separator = ','
    elif ligand_file_format == 'tsv':
        field_separator = '\t'
    else:
        field_separator = None

    if ligand_file_format == 'pkl':
        df = pd.read_pickle(ligand_list_file)
    else:
        df = pd.read_csv(ligand_list_file, sep=field_separator)

#    num_data = df.shape[0]
    fkey = df.keys()[0]
    if fkey.startswith('#'):
        df.rename(columns={fkey: fkey[1:]}, inplace=True)
    smiles_list = df[['MOL_ID', 'SMILES']].values.tolist()

#    smiles_list = smiles_list[0:10]
    docking_vina = DockingVina(docking_params)

    result_dict = docking_vina.predict(smiles_list)
    docking_score_list = result_dict['docking']
    docking_min = [x[0] for x in docking_score_list]
    df['Docking1'] = docking_min
    df['Docking'] = docking_score_list
    if docking_params['check_metal_bind']:
        num_metal_bind_atom_list = result_dict['num_metal_bind_atom']
        df['Metal_bind'] = num_metal_bind_atom_list
    if docking_params['rescoring']:
        rescoring = result_dict['docking_re']
        df['Docking_re'] = rescoring

    use_my_module = docking_params['use_my_module']
    my_module_path = docking_params['my_module_path']

    if use_my_module:
        my_module_dir = os.path.dirname(my_module_path)
        sys.path.append(my_module_dir)
        import my_module
        my_module.my_score_to_df(df, docking_params, result_dict)

    sep = field_separator
    if sep == '\s+':
        sep = ' '

    if output_file.strip().split('.')[-1] == 'pkl':
        df.to_pickle(output_file)
    else:
        df.to_csv(output_file, sep=sep, float_format='%.3f', index=False)


if __name__ == "__main__":
    main()
