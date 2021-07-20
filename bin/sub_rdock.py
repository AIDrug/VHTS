#!/usr/bin/env python
import sys
import os
import argparse
import vhts.pyrdock as pyrdock
from filelock import FileLock
import pandas as pd


def get_job_from_list(list_dir):
    list_file = list_dir + '/list.txt'
    if not os.path.exists(list_file):
        job_idx = None
        return job_idx
    freeze_lock = FileLock('{}.lock'.format(list_file))
    with freeze_lock.acquire(timeout=30):
        with open(list_file, 'r') as fp:
            lines = fp.readlines()
        if len(lines) == 0:
            job_idx = None
        else:
            job_idx = lines[0].strip()
            with open(list_file, 'w') as fp:
                for line in lines[1:]:
                    fp.write(line)
    return job_idx


def set_job_from_list(job_idx, list_dir):
    list_file = list_dir + '/list.txt'
    freeze_lock = FileLock('{}.lock'.format(list_file))
    with freeze_lock.acquire(timeout=30):
        if os.path.exists(list_file):
            with open(list_file, 'r') as fp:
                lines = fp.readlines()
        else:
            lines = list()
        with open(list_file, 'w') as fp:
            for line in lines:
                fp.write(line)
            line = job_idx + '\n'
            fp.write(line)
    return


def remove_job_from_list(job_idx, list_dir):
    list_file = list_dir + '/list.txt'
    freeze_lock = FileLock('{}.lock'.format(list_file))
    with freeze_lock.acquire(timeout=30):
        if os.path.exists(list_file):
            with open(list_file, 'r') as fp:
                lines = fp.readlines()
        with open(list_file, 'w') as fp:
            for line in lines:
                if line.strip() != job_idx:
                    fp.write(line)
    return


def get_and_set_job(params_dict):
    todo_dir = params_dict['todo_dir']
    current_dir = params_dict['current_dir']
    job_idx = get_job_from_list(todo_dir)
    if job_idx is None:
        return job_idx
    set_job_from_list(job_idx, current_dir)
    smi_file_format = params_dict['smi_file_format']
    job_todo_file = todo_dir + '/' + job_idx + '.' + smi_file_format
    job_current_file = current_dir + '/' + job_idx + '.' + smi_file_format

    os.replace(job_todo_file, job_current_file)
    return job_idx


def run_docking(job_idx, docking, params_dict):
    current_dir = params_dict['current_dir']
    done_dir = params_dict['done_dir']
    field_separator = params_dict['field_separator']
    smi_file_format = params_dict['smi_file_format']
    job_current_file = current_dir + '/' + job_idx + '.' + smi_file_format
    job_done_file = done_dir + '/' + job_idx + '.' + smi_file_format
    if smi_file_format == 'pkl':
        df = pd.read_pickle(job_current_file)
    else:
        df = pd.read_csv(job_current_file, sep=field_separator, header=0)

#    num_data = df.shape[0]
#    df.rename(columns={0: 'MOL_IDX', 1: 'MOL_ID', 2: 'SMILES'}, inplace=True)
    smiles_list = df[['MOL_ID', 'SMILES']].values.tolist()
    result_dict = docking.predict(smiles_list)
    affinity_list = result_dict['docking']
    docking_min = [x[0] for x in affinity_list]
#    docking = [x for x in affinity_list]
    docking = affinity_list
    df['Docking1'] = docking_min
    df['Docking'] = docking
    if params_dict['rescoring']:
        rescoring = result_dict['docking_re']
        df['Docking_re'] = rescoring

    use_my_module = params_dict['use_my_module']
    my_module_path = params_dict['my_module_path']
    docking_params = params_dict['docking_params']

    if use_my_module:
        my_module_dir = os.path.dirname(my_module_path)
        sys.path.append(my_module_dir)
        import my_module
        my_module.my_score_to_df(df, docking_params, result_dict)

    sep = field_separator
    if sep == '\s+':
        sep = ' '

    if smi_file_format == 'pkl':
        df.to_pickle(job_done_file)
    else:
        df.to_csv(job_done_file, sep=sep, float_format='%.3f',
                  header=True, index=False)

    return


def move_done(job_idx, params_dict):
    current_dir = params_dict['current_dir']
    done_dir = params_dict['done_dir']
    remove_job_from_list(job_idx, current_dir)
    set_job_from_list(job_idx, done_dir)
    smi_file_format = params_dict['smi_file_format']
    job_current_file = current_dir + '/' + job_idx + '.' + smi_file_format
    os.remove(job_current_file)
    return job_idx


def working(docking, params_dict):
    pid = os.getpid()
    out_log = params_dict['out_log']
    log_file = params_dict['log_file']

    line_out = 'Start sub_dock pid: %d' % (pid)
    if out_log == 'file':
        fp_log = open(log_file, 'w')
        fp_log.write(line_out + '\n')
        fp_log.flush()
    elif out_log == 'print':
        print(line_out, flush=True)

    while True:
        job_idx = get_and_set_job(params_dict)
        line_out = 'get a job: %s' % job_idx
        if out_log == 'file':
            fp_log.write(line_out + '\n')
            fp_log.flush()

        elif out_log == 'print':
            print(line_out, flush=True)
        if job_idx is None:
            line_out = 'End sub_dock pid %d' % pid
            if out_log == 'file':
                fp_log.write(line_out + '\n')
                fp_log.flush()

            elif out_log == 'print':
                print(line_out, flush=True)
            break
        run_docking(job_idx, docking, params_dict)
        move_done(job_idx, params_dict)
        line_out = 'done job: %s' % job_idx
        if out_log == 'file':
            fp_log.write(line_out + '\n')
            fp_log.flush()

        elif out_log == 'print':
            print(line_out, flush=True)
    if out_log == 'file':
        fp_log.close()

    return


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

    parser = argparse.ArgumentParser(description='worker for docking')
    parser.add_argument('--work_dir', type=str, required=False,
                        default='workflow', help='workflow directory')
    parser.add_argument('-s', '--smi_file_format', type=str, required=False,
                        default='pkl', help='pkl (default), txt, csv, tsv')
    parser.add_argument('--out_log', type=str, required=False,
                        default=None,
                        help='print : screen, or file : sub_dock_$PID.log' +
                        'default: do not print output')

    parser_arg(parser)

    args, docking_params = arg_to_params(parser)
    if len(sys.argv) < 2:
        parser.print_usage()
        sys.exit()
    if args.dock_config is None:
        parser.print_usage()
        print('dock_config is missing')
        sys.exit()

    work_dir = args.work_dir
    todo_dir = work_dir + '/todo'
    current_dir = work_dir + '/current'
    done_dir = work_dir + '/done'
    smi_file_format = args.smi_file_format
    if smi_file_format == 'txt':
        field_separator = '\s+'
    elif smi_file_format == 'csv':
        field_separator = ','
    elif smi_file_format == 'tsv':
        field_separator = '\t'
    else:
        field_separator = None

    out_log = args.out_log
    pid = os.getpid()
    log_file = 'sub_dock_%d.log' % (pid)

    docking = pyrdock.Docking(docking_params)

    params_dict = dict()
    params_dict['work_dir'] = work_dir
    params_dict['todo_dir'] = todo_dir
    params_dict['current_dir'] = current_dir
    params_dict['done_dir'] = done_dir
    params_dict['field_separator'] = field_separator
    params_dict['smi_file_format'] = smi_file_format
    params_dict['out_log'] = out_log
    params_dict['log_file'] = log_file

    params_dict['rescoring'] = docking_params['rescoring']
    params_dict['use_my_module'] = docking_params['use_my_module']
    params_dict['my_module_path'] = docking_params['my_module_path']
    params_dict['docking_params'] = docking_params

    working(docking, params_dict)


if __name__ == "__main__":
    main()
