#!/usr/bin/env python
import sys
import os
import numpy as np
from multiprocessing import Manager
from multiprocessing import Process
from multiprocessing import Queue
import subprocess
import argparse
import pandas as pd
from pdbtools import ligand_tools
from vhts import sdtether


class Docking(object):
    """
    python module for Docking
    """

    def __init__(self, docking_params):
        """
            Construction Docking object
        """
        self.docking_program = docking_params['docking_program']
        dp_tmp = self.docking_program.lower()
        if (dp_tmp.find('vina') != -1) or (dp_tmp.find('smina') != -1):
            self.dp = 'vina'
        elif (dp_tmp.find('rbdock') != -1):
            self.dp = 'rdock'
        else:
            self.dp = 'etc'
        if 'exhaustiveness' not in docking_params:
            self.exhaustiveness = 10
        else:
            self.exhaustiveness = docking_params['exhaustiveness']
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

        self.rescoring = docking_params['rescoring']
        self.rescoring_program = docking_params['rescoring_program']
        self.rescoring_config_file = docking_params['rescoring_config_file']

        self.tether_SMARTS = docking_params['tether_SMARTS']
        self.tether_ref_lig = docking_params['tether_ref_lig']
        self.tether_ref_coor_file = docking_params['tether_ref_coor_file']

        self.tether_docking = False
        if self.tether_SMARTS is not None:
            if self.dp != 'rdock':
                print('Tether scaffold docking is only supported by rdock.')
                sys.exit()
            self.smarts = sdtether.gen_smarts(self.tether_SMARTS)
            self.ref_match_coords = list()
            if self.tether_ref_coor_file is not None:
                self.ref_match_coords = sdtether.read_tether_coords(
                    self.tether_ref_coor_file)
            elif self.tether_ref_lig is not None:
                self.ref_match_coords = sdtether.find_tether_ref_coords(
                    self.tether_ref_ligand, self.tether_SMARTS)

            if len(self.ref_match_coords) >= 1:
                self.tether_docking = True
            else:
                print('No reference matching coordinates.')
                sys.exit()

    def docking_vina(self, ligand_file, docking_pdbqt_file, docking_log_file):
        """
            run_docking program using subprocess
            input :
                ligand_file
                docking_pdbqt_file
            output :
                affinity list for a input molecule
        """

        run_line = '%s' % self.docking_program
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

    def docking_vina_score_only(self, ligand_file):
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

    def docking_rdock(self, ligand_file, docking_file, docking_log_file):
        """
            run_docking program using subprocess
            input :
                ligand_file
                docking_pdbqt_file
            output :
                affinity list for a input molecule
        """

        docking_prefix = '.'.join(docking_file.strip().split('.')[:-1])
        run_line = '%s' % self.docking_program
        run_line += ' -r %s' % self.dock_config_file
        run_line += ' -p dock.prm'
        run_line += ' -n %d' % self.exhaustiveness
        run_line += ' -i %s' % ligand_file
        run_line += ' -o %s' % docking_prefix

#        run_line2 = 'sdsort -n -fSCORE %s.sd' % (docking_prefix)
        run_line2 = 'sdsort -n -fSCORE.INTER %s.sd' % (docking_prefix)

        e = None
        try:
            result = subprocess.check_output(run_line.split(),
                                             stderr=subprocess.STDOUT,
                                             timeout=self.timeout_dock,
                                             universal_newlines=True)
            if self.output_save:
                fp = open(docking_log_file, 'w')
                fp.write(result)
                fp.close()

            result2 = subprocess.check_output(run_line2.split(),
                                              universal_newlines=True)
            fp = open(docking_file, 'w')
            fp.write(result2)
            fp.close()

        except Exception as e:
            return [99.999], e

        affinity_list = list()
        out_lines = result2.split('\n')
        check_score = False
        for line in out_lines:
            if line[0:16] == '>  <SCORE.INTER>':
#            if line[0:10] == '>  <SCORE>':
                check_score = True
                continue
            if check_score is True:
                affinity = float(line)
                affinity_list += [affinity]
                check_score = False
                continue
        if len(affinity_list) == 0:
            e = 'WARNING: Could not find any conformations.'
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
            smi_p = ligand_tools.ligand_preparation(smi, self.neutralize,
                                                    self.pH)
        else:
            smi_p = smi
        if not self.output_save:
            ligand_sdf_file = '%s/ligand_%d.sdf' % (self.gen3d_dir, pid)
            ligand_pdb_file = '%s/ligand_%d.pdb' % (self.gen3d_dir, pid)
            ligand_pdbqt_file = '%s/ligand_%s.pdbqt' % (self.gen3d_dir, pid)
            docking_pdbqt_file = '%s/dock_%d.pdbqt' % (
                self.dock_dir, pid)
            docking_log_file = '%s/dock_%d.log' % (self.dock_dir, pid)
            docking_pdb_file = '%s/dock_%s.pdb' % (self.dock_dir, pid)
            docking_sdf_file = '%s/dock_%s.sdf' % (self.dock_dir, pid)

            out_dock_dir1 = None
        else:
            mol_id2 = mol_id[0:self.tlen]
            out_gen3d_dir1 = self.gen3d_dir + "/" + mol_id2
            if not os.path.exists(out_gen3d_dir1):
                try:
                    os.makedirs(out_gen3d_dir1)
                except FileExistsError as e:
                    print(e, flush=True)
            ligand_sdf_file = '%s/ligand_%s.sdf' % (out_gen3d_dir1, mol_id)
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
            docking_sdf_file = '%s/dock_%s.sdf' % (out_dock_dir1, mol_id)

        if self.dp == 'vina':
            e = ligand_tools.gen_3d(smi_p, ligand_pdb_file, mol_id=mol_id,
                                    timeout=self.timeout_gen3d)
            if e is not None:
                e2 = ligand_tools.gen_3d(smi_p, ligand_pdb_file, mol_id=mol_id,
                                         timeout=self.timeout_gen3d)
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

        elif self.dp == 'rdock':
            e = ligand_tools.gen_3d(smi_p, ligand_sdf_file, mol_id=mol_id,
                                    timeout=self.timeout_gen3d)
            if e is not None:
                e2 = ligand_tools.gen_3d(smi_p, ligand_sdf_file, mol_id=mol_id,
                                         timeout=self.timeout_gen3d)
                if e2 is not None:
                    print(e2, 'gen_3d', idx, mol_id, smi_p, flush=True)
                    docking_score = np.array([99.999], dtype=np.float32)
                    result_dict['docking'] = docking_score
                    return result_dict

            if self.tether_docking:
                e = sdtether.tether_ligand(ligand_sdf_file, self.smarts,
                                           self.ref_match_coords)
                if e is not None:
                    docking_score = np.array([99.999], dtype=np.float32)
                    result_dict['docking'] = docking_score
                    return result_dict

        if self.dp == 'vina':
            docking_score, e = self.docking_vina(ligand_pdbqt_file,
                                                 docking_pdbqt_file,
                                                 docking_log_file)
        elif self.dp == 'rdock':
            docking_score, e = self.docking_rdock(ligand_sdf_file,
                                                  docking_sdf_file,
                                                  docking_log_file)

        docking_score = np.array(docking_score, dtype=np.float32)
        if e is not None:
            docking_score = [99.999]
            result_dict['docking'] = docking_score
            print(e, 'docking', idx, mol_id, smi_p, flush=True)
            return result_dict
        result_dict['docking'] = docking_score
        if self.output_save or self.rescoring or self.use_my_module:
            if self.dp == 'vina':
                ligand_tools.pdbqt_to_pdb_ref(docking_pdbqt_file,
                                              docking_pdb_file,
                                              ligand_pdb_file)
            elif self.dp == 'rdock':
                ligand_tools.obabel_rewrite(docking_sdf_file,
                                            docking_pdb_file, option=' -h')

        if self.rescoring:
            docking_rescore, e = self.docking_vina_score_only(docking_pdb_file)
            docking_rescore = np.array(docking_rescore, dtype=np.float32)
            if e is not None:
                docking_rescore = np.array([99.999], dtype=np.float32)
                print(e, 're-scoring', idx, mol_id, smi_p, flush=True)
            result_dict['docking_re'] = docking_rescore

        if self.use_my_module:
            self.my_class.simulation_process(self, idx, mol_id, smi, smi_p,
                                             pid, out_dock_dir1,
                                             docking_pdb_file, result_dict)

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
        if self.rescoring:
            docking_re_list = list()

        for key in range(num_data):
            if key in keys:
                result_dict0 = return_dict[key]
                if 'docking' in result_dict0:
                    docking_score = result_dict0['docking']
                else:
                    docking_score = np.array([99.999], dtype=np.float32)

                if self.rescoring:
                    if 'docking_re' in result_dict0:
                        docking_re = result_dict0['docking_re']
                    else:
                        docking_re = np.array([99.999], dtype=np.float32)

            else:
                docking_score = np.array([99.999], dtype=np.float32)
                if self.rescoring:
                    docking_re = np.array([99.999], dtype=np.float32)

            docking_score_list += [docking_score]
            if self.rescoring:
                docking_re_list += [docking_re]

        result_dict['docking'] = docking_score_list
        if self.rescoring:
            result_dict['docking_re'] = docking_re_list

        if self.use_my_module:
            self.my_class.predict(self, smiles_list, result_dict, return_dict)

        return result_dict


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
    # docking parameter

    parser.register('action', 'extend', ExtendAction)
    parser.add_argument('--arg_file', type=open, required=False, default=None,
                        action=LoadFromConfig, help='argment file')
    parser.add_argument('--dock_config', type=str, required=False,
                        default=None, help='docking config file ')
    parser.add_argument('-v', '--docking_program', type=str, required=False,
                        default='rbdock',
                        help='select rdock, rbdock')
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
    parser.add_argument('--rescoring_program', type=str, required=False,
                        default='smina', help='smina path')
    parser.add_argument('--rescoring_config', type=str, required=False,
                        default=None, help='docking config file for rescoring')

    parser.add_argument('--tether_ref_lig', type=str, required=False,
                        default=None,
                        help='reference ligand for tether docking')
    parser.add_argument('--tether_SMARTS', type=str, required=False,
                        default=None, help='SMARTS pattern for tether docking')
    parser.add_argument('--tether_ref_coor_file', type=str, required=False,
                        default=None,
                        help='reference coordinate file for tether docking')

    parser.add_argument('--exhaustiveness', type=int, required=False,
                        default=10,
                        help='exhaustiveness for rdock')

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

    docking_program = args.docking_program
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

    tether_ref_lig = args.tether_ref_lig
    tether_SMARTS = args.tether_SMARTS
    tether_ref_coor_file = args.tether_ref_coor_file

    exhaustiveness = args.exhaustiveness

    rescoring = False
    rescoring_config_file = args.rescoring_config
    rescoring_program = args.rescoring_program
    if rescoring_config_file is not None:
        rescoring = True

    docking_params = dict()
    docking_params['docking_program'] = docking_program
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
    docking_params['rescoring'] = rescoring
    docking_params['rescoring_program'] = rescoring_program
    docking_params['rescoring_config_file'] = rescoring_config_file
    docking_params['tether_ref_lig'] = tether_ref_lig
    docking_params['tether_SMARTS'] = tether_SMARTS
    docking_params['tether_ref_coor_file'] = tether_ref_coor_file
    docking_params['exhaustiveness'] = exhaustiveness

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
    docking = Docking(docking_params)

    result_dict = docking.predict(smiles_list)
    docking_score_list = result_dict['docking']
    docking_min = [x[0] for x in docking_score_list]
    df['Docking1'] = docking_min
    df['Docking'] = docking_score_list
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
