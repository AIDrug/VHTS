import sys
import numpy as np
from vhts import metal_bind


class my_class(object):
    def __init__(self, docking_params):
        self.check_metal_bind = docking_params['check_metal_bind']
        if self.check_metal_bind:
            self.metal_coor = docking_params['metal_coor']
            self.metal_bind_cutoff = docking_params['metal_bind_cutoff']

    def simulation_process(self, idx, mol_id, smi, pid, smi_p,
                           out_dock_dir1, docking_pdb_file, result_dict):

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

    def predict(self, smiles_list, result_dict, return_dict):
        data = list(enumerate(smiles_list))
        num_data = len(data)

        keys = sorted(return_dict.keys())

        if self.check_metal_bind:
            num_metal_bind_atom_list = list()

        for key in range(num_data):
            if key in keys:
                result_dict0 = return_dict[key]

                if self.check_metal_bind:
                    if 'num_metal_bind_atom' in result_dict0:
                        num_metal_bind_atom = result_dict0['num_metal_bind_atom']
                    else:
                        num_metal_bind_atom = np.array([0], dtype=np.float32)
            else:
                if self.check_metal_bind:
                    num_metal_bind_atom = np.array([0], dtype=np.float32)
            if self.check_metal_bind:
                num_metal_bind_atom_list += [num_metal_bind_atom]
        if self.check_metal_bind:
            result_dict['num_metal_bind_atom'] = num_metal_bind_atom_list
        return result_dict


def parser_arg(parser):
    parser.add_argument('--metal_coor', type=str, default=None, required=False,
                        help='position of metal_ion,' +
                        ' example: --metal_coor="1.0,-1.0,0.0" default: None')
    parser.add_argument('--metal_cutoff', type=float, required=False,
                        default=3.0,
                        help='metal ion - HDA cutoff distance,')


def arg_to_params(parser, docking_params):

    args = parser.parse_args()

    check_metal_bind = False
    if args.metal_coor is not None:
        check_metal_bind = True

        metal_coor = np.array(args.metal_coor.strip(
            '"').split(','), dtype=np.float32)
        if metal_coor.shape[0] != 3:
            print('metal coordinate is strange', args.metal_coor)
            sys.exit()
        metal_bind_cutoff = args.metal_cutoff

        docking_params['check_metal_bind'] = check_metal_bind
        docking_params['metal_coor'] = metal_coor
        docking_params['metal_bind_cutoff'] = metal_bind_cutoff
    return docking_params


def my_score_to_df(df, docking_params, result_dict):

    if docking_params['check_metal_bind']:
        num_metal_bind_atom_list = result_dict['num_metal_bind_atom']
        df['Metal_bind'] = num_metal_bind_atom_list
