import numpy as np
from pifinder import pifinder


class my_class(object):
    def __init__(self, docking_params):
        self.use_piscore = docking_params['use_piscore']
        if self.use_piscore:
            self.pharma = docking_params['pharma']
            self.use_pinfo = self.pharma.use_pinfo
            self.vpiscore_w1 = docking_params['vpiscore_w1']
            self.vpiscore_w2 = docking_params['vpiscore_w2']

    def simulation_process(self, idx, mol_id, smi, pid, smi_p,
                           out_dock_dir1, docking_pdb_file, result_dict):
        if self.use_piscore:
            if self.output_save:
                pf_ligand_file = '%s/pf_%s.txt' % (out_dock_dir1, mol_id)
                interaction_file = '%s/pinter_%s.txt' % (out_dock_dir1, mol_id)
            else:
                pf_ligand_file = None
                interaction_file = None
            try:
                result_piscore = self.pharma.cal_piscore(
                    docking_pdb_file, pf_ligand_file, interaction_file)
                total_count_dict, total_count_info_dict = result_piscore
                keys = list(total_count_dict.keys())
                piscore = list()
                if self.use_pinfo:
                    pinfo = list()
                for key in keys:
                    count_dict = total_count_dict[key]
                    pis = count_dict['Score']
                    piscore += [pis]
                    if self.use_pinfo:
                        count_info_dict = total_count_info_dict[key]
                        pin = count_info_dict['Score']
                        pinfo += [pin]
                piscore = np.array(piscore, dtype=np.float32)
                if self.use_pinfo:
                    pinfo = np.array(pinfo, dtype=np.float32)

            except Exception as e:
                print(e, 'piscore', idx, mol_id, smi_p, flush=True)
                piscore = np.array([0], dtype=np.float32)
                if self.use_pinfo:
                    pinfo = np.array([0], dtype=np.float32)

            result_dict['piscore'] = piscore
            if self.rescoring:
                docking_rescore = result_dict['docking_re']
                vpiscore = (self.vpiscore_w1 * (-docking_rescore)
                            + (1.0 - self.vpiscore_w1) * piscore)
            else:
                docking_score = result_dict['docking']
                vpiscore = (self.vpiscore_w1 * (-docking_score)
                            + (1.0 - self.vpiscore_w1) * piscore)

            result_dict['vpiscore'] = vpiscore

            if self.use_pinfo:
                result_dict['pinfo'] = pinfo
                if self.rescoring:
                    docking_rescore = result_dict['docking_re']
                    vpiscore_info = (self.vpiscore_w1 * (-docking_rescore)
                                     + (1.0 - self.vpiscore_w1)
                                     * (self.vpiscore_w2 * piscore
                                        + (1.0 - self.vpiscore_w1) * pinfo))
                else:
                    docking_score = result_dict['docking']
                    vpiscore_info = (self.vpiscore_w1 * (-docking_score)
                                     + (1.0 - self.vpiscore_w1)
                                     * (self.vpiscore_w2 * piscore
                                        + (1.0 - self.vpiscore_w1) * pinfo))

                result_dict['vpiscore_info'] = vpiscore_info

    def predict(self, smiles_list, result_dict, return_dict):
        data = list(enumerate(smiles_list))
        num_data = len(data)

        keys = sorted(return_dict.keys())
        if self.use_piscore:
            piscore_list = list()
            vpiscore_list = list()
            if self.use_pinfo:
                pinfo_list = list()
                vpiscore_info_list = list()

        for key in range(num_data):
            if key in keys:
                result_dict0 = return_dict[key]
                if self.use_piscore:
                    if 'piscore' in result_dict0:
                        piscore = result_dict0['piscore']
                    else:
                        piscore = np.array([0], dtype=np.float32)

                    if 'vpiscore' in result_dict0:
                        vpiscore = result_dict0['vpiscore']
                    else:
                        vpiscore = np.array([0], dtype=np.float32)
                    if self.use_pinfo:
                        if 'pinfo' in result_dict0:
                            pinfo = result_dict0['pinfo']
                        else:
                            pinfo = np.array([0], dtype=np.float32)
                        if 'vpiscore_info' in result_dict0:
                            vpiscore_info = result_dict0['vpiscore_info']
                        else:
                            vpiscore_info = np.array([0], dtype=np.float32)
            else:
                if self.use_piscore:
                    piscore = np.array([0], dtype=np.float32)
                    vpiscore = np.array([0], dtype=np.float32)
                    if self.use_pinfo:
                        pinfo = np.array([0], dtype=np.float32)
                        vpiscore_info = np.array([0], dtype=np.float32)
            if self.use_piscore:
                piscore_list += [piscore]
                vpiscore_list += [vpiscore]
                if self.use_pinfo:
                    pinfo_list += [pinfo]
                    vpiscore_info_list += [vpiscore_info]
        if self.use_piscore:
            result_dict['piscore_list'] = piscore_list
            result_dict['vpiscore_list'] = vpiscore_list
            if self.use_pinfo:
                result_dict['pinfo_list'] = pinfo_list
                result_dict['vpiscore_info_list'] = vpiscore_info_list
        return result_dict


def parser_arg(parser):
    parser.add_argument('--piscore_receptor', required=False,
                        default=None,
                        help='input receptor pdb file for piscore')
    parser.add_argument('--pf_receptor', required=False,
                        default=None,
                        help='output for pharmacophoric feature of receptor')
    parser.add_argument('--pinfo_ligand', required=False,
                        default=None, help='template ligand file for pinfo')
    parser.add_argument('--pf_receptor_info', required=False, default=None,
                        help='output for template feature of receptor')
    parser.add_argument('--pi_cutoff', type=float, default=6.5, required=False,
                        help='pharmacophoric interaction cutoff distance')
    parser.add_argument('--vpiscore_w1', type=float, required=False,
                        default=0.2,
                        help='weight_1 of VPIscore\n' +
                        'VPIscore = w1*vina_score + (1-w1)*piscore')
    parser.add_argument('--vpiscore_w2', type=float, required=False,
                        default=0.8,
                        help='weight_2 of VPIscore_info\n' +
                        'VPIscore_info = w1*vina_score + (1-w1)*(w2*piscore'
                        + '(1-w2)*pinfo)')
    parser.add_argument('--include_hydrophobic', action='store_true',
                        required=False,
                        help='include hydrophobic feature for template')


def arg_to_params(parser, docking_params):

    args = parser.parse_args()
    use_piscore = False
    if args.pinfo_ligand is not None or args.pf_receptor_info is not None:
        use_piscore = True
        pharma = pifinder.set_pifinder(args)
        docking_params['use_piscore'] = use_piscore
        docking_params['pharma'] = pharma
        docking_params['use_pinfo'] = pharma.use_pinfo
        docking_params['vpiscore_w1'] = args.vpiscore_w1
        docking_params['vpiscore_w2'] = args.vpiscore_w2
    return docking_params


def my_score_to_df(df, docking_params, result_dict):
    if docking_params['use_piscore']:
        piscore = result_dict['piscore_list']
        vpiscore = result_dict['vpiscore_list']
        vpiscore_1 = [x.max() for x in vpiscore]
        df['VPIscore1'] = vpiscore_1

        if docking_params['use_pinfo']:
            pinfo = result_dict['pinfo_list']
            vpiscore_info = result_dict['vpiscore_info_list']
            vpiscore_info_1 = [x.max() for x in vpiscore_info]

            df['VPIscore_info1'] = vpiscore_info_1

        df['VPIscore'] = vpiscore
        if docking_params['use_pinfo']:
            df['VPIscore_info'] = vpiscore_info
        df['PIscore'] = piscore
        if docking_params['use_pinfo']:
            df['Pinfo'] = pinfo
