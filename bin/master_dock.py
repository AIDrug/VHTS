#!/usr/bin/env python
import sys
import os
import numpy as np
import argparse
import time
from filelock import FileLock
import pandas as pd
import subprocess


def get_job_list(list_dir):
    job_idx_list = list()
    list_file = list_dir + '/list.txt'
    if not os.path.exists(list_file):
        return job_idx_list
    freeze_lock = FileLock('{}.lock'.format(list_file))
    with freeze_lock.acquire(timeout=30):
        with open(list_file, 'r') as fp:
            lines = fp.readlines()
        for line in lines:
            job_idx_list += [line.strip()]
    job_idx_list2 = list()
    for job_idx in job_idx_list:
        x = job_idx.split('_')
        job_idx_list2 += [[int(x[0]), int(x[1])]]
    job_idx_list2 = sorted(job_idx_list2)
    job_idx_list3 = ['%d_%d' % (x[0], x[1]) for x in job_idx_list2]

    return job_idx_list3


def set_job_list(job_idx_list, list_dir):
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
            for job_idx in job_idx_list:
                line = job_idx + '\n'
                fp.write(line)

    return


def data_sampling(params_dict, iteration):

    field_separator = params_dict['field_separator']
    smi_file_format = params_dict['smi_file_format']
    sampling_method = params_dict['sampling_method']
    num_sample = params_dict['num_sample']
    file_size = params_dict['file_size']
    master_dir = params_dict['master_dir']
    todo_dir = params_dict['todo_dir']

    remain_smi_file = master_dir + '/' + 'remain.' + smi_file_format

    try:
        if smi_file_format == 'pkl':
            df_remain = pd.read_pickle(remain_smi_file)
        else:
            df_remain = pd.read_csv(
                remain_smi_file, sep=field_separator, header=0)

    except pd.errors.EmptyDataError:
        return False
#    0: 'MOL_IDX', 1: 'MOL_ID', 2: 'SMILES' 3:, 'Pred', 4: 'uncertainty'
    fkey = df_remain.keys()[0]
    if fkey.startswith('#'):
        df_remain.rename(columns={fkey: fkey[1:]}, inplace=True)

    num_remain = df_remain.shape[0]
    if num_remain == 0:
        return False, False
    if sampling_method == 'sequential':
        df_sample = df_remain[0:num_sample]
        df_remain_new = df_remain[num_sample:]
    elif sampling_method == 'random' or iteration == 0:
        index = df_remain.index.values
        np.random.shuffle(index)
        idx_sam = np.sort(index[0:num_sample])
        idx_rem = np.sort(index[num_sample])
        df_sample = df_remain.loc[idx_sam]
        df_remain_new = df_remain.loc[idx_rem]
    elif sampling_method == 'score':
        index = df_remain['Pred'].sort_values().index
        idx_sam = np.sort(index[0:num_sample])
        idx_rem = np.sort(index[num_sample])
        df_sample = df_remain.loc[idx_sam]
        df_remain_new = df_remain.loc[idx_rem]
    elif sampling_method == 'uncertainty':
        index = df_remain['Uncertainty'].sort_values(ascending=False).index
        idx_sam = np.sort(index[0:num_sample])
        idx_rem = np.sort(index[num_sample])
        df_sample = df_remain.loc[idx_sam]
        df_remain_new = df_remain.loc[idx_rem]

    sample_smi_file = '%s/sample_%d.%s' % (master_dir, iteration,
                                           smi_file_format)
    sep = field_separator
    if sep == '\s+':
        sep = ' '

    if smi_file_format == 'pkl':
        df_sample.to_pickle(sample_smi_file)
        df_remain_new.to_pickle(remain_smi_file)

    else:
        df_sample.to_csv(sample_smi_file, sep=sep, float_format='%.3f',
                         header=True, index=False)
        df_remain_new.to_csv(remain_smi_file, sep=sep,
                             float_format='%.3f', index=False)

    todo_list = list()

    num_data = df_sample.shape[0]
    num_file = int(np.ceil(num_data/file_size))
    for idx in range(0, num_file):
        job_idx = '%d_%d' % (iteration, idx)
        ini = idx * file_size
        fin = (idx+1) * file_size
        df_todo = df_sample[ini:fin]
        job_todo_file = todo_dir + '/' + job_idx + '.' + smi_file_format

        if smi_file_format == 'pkl':
            df_todo.to_pickle(job_todo_file)
        else:
            df_todo.to_csv(job_todo_file, sep=sep, float_format='%.3f',
                           header=True, index=False)
        todo_list += [job_idx]

    set_job_list(todo_list, master_dir)
    set_job_list(todo_list, todo_dir)
    num_remain_new = df_remain_new.shape[0]
    check_remain = True
    if num_remain_new == 0:
        check_remain = False

    return True, check_remain


def write_stage(stage_list, stage_file):
    fp_stage = open(stage_file, 'w')
    for stage in stage_list:
        iteration, state = stage
        line_out = '%d %s\n' % (iteration, state)
        fp_stage.write(line_out)
    fp_stage.close()


def check_done(params_dict, iteration):

    check = False
    master_dir = params_dict['master_dir']
    done_dir = params_dict['done_dir']
#    mlsync = params_dict['mlsync']

    master_job_idx_list = get_job_list(master_dir)
    done_job_idx_list = get_job_list(done_dir)
    master_job_i = list()
    done_job_i = list()

    for job_idx in master_job_idx_list:
        lis = job_idx.strip().split('_')
        idx = int(lis[1])
        iteration0 = int(lis[0])
        if iteration == iteration0:
            master_job_i += [idx]
    for job_idx in done_job_idx_list:
        lis = job_idx.strip().split('_')
        idx = int(lis[1])
        iteration0 = int(lis[0])
        if iteration == iteration0:
            done_job_i += [idx]

    master_job = set(master_job_i)
    done_job = set(done_job_i)
    running_job = master_job - done_job
    if len(running_job) == 0:
        check = True
    return check


def gather_result(params_dict, iteration):

    master_dir = params_dict['master_dir']
    done_dir = params_dict['done_dir']
    field_separator = params_dict['field_separator']
    smi_file_format = params_dict['smi_file_format']

#    head_dict = {0: 'MOL_IDX', 1: 'MOL_ID',
#                 2: 'SMILES', 3: 'Docking1', 4: 'Docking'}
    docking_smi_file = master_dir + '/' + 'docking.' + smi_file_format

    all_data = list()
    if os.path.exists(docking_smi_file):
        if smi_file_format == 'pkl':
            df_docking = pd.read_pickle(docking_smi_file)
        else:
            df_docking = pd.read_csv(
                docking_smi_file, sep=field_separator, header=0)
        all_data.append(df_docking)
    done_job_idx_list = get_job_list(done_dir)
    for job_idx in done_job_idx_list:
        lis = job_idx.strip().split('_')
        iteration0 = int(lis[0])
        if iteration != iteration0:
            continue
        job_done_file = done_dir + '/' + job_idx + '.' + smi_file_format
        if smi_file_format == 'pkl':
            df_done = pd.read_pickle(job_done_file)
        else:
            df_done = pd.read_csv(job_done_file, sep=field_separator, header=0)
#        df_done.rename(columns=head_dict, inplace=True)
        all_data.append(df_done)
    df = pd.concat(all_data, axis=0, ignore_index=True)
    sep = field_separator
    if sep == '\s+':
        sep = ' '
    if smi_file_format == 'pkl':
        df.to_pickle(docking_smi_file)
    else:
        df.to_csv(docking_smi_file, sep=sep, float_format='%.3f', index=False)
    return


def gapjil(params_dict):
    '''
        master gapjil
    '''
    sleep_cycle = params_dict['sleep_cycle']
#    sampling_method = params_dict['sampling_method']
    mlsync = params_dict['mlsync']
    master_dir = params_dict['master_dir']
    auto_sub = params_dict['auto_sub']
    sub_dock_script = params_dict['sub_dock_script']
    num_sub = params_dict['num_sub']

    stage_file = master_dir + '/' + 'stage.txt'
    stage_list = list()
    iteration = 0

    check_docking = True
    check_remain = True
    if os.path.exists(stage_file):
        stage_list = [x.strip().split() for x in open(stage_file)]
        if len(stage_list) != 0:
            iteration = int(stage_list[-1][0])
            if stage_list[-1][1] != 'done':
                check_docking = False
        line_out = 'restart iteration: %d docking_done: %s' % (
            iteration, check_docking)
        print(line_out, flush=True)
    while True:
        if not check_docking:
            check_docking = check_done(params_dict, iteration)
            if check_docking:
                stage_list[iteration] = [iteration, 'done']
                write_stage(stage_list, stage_file)
                gather_result(params_dict, iteration)
                line_out = 'end iteration: %d docking_done: %s' % (
                    iteration, check_docking)
                print(line_out, flush=True)
                iteration += 1
                # do ml
        if check_docking and not check_remain:
            line_out = 'End simulation'
            print(line_out, flush=True)
            break
        if check_docking and check_remain:
            check_sampling, check_remain = data_sampling(
                params_dict, iteration)
            if not check_sampling:
                line_out = 'End simulation'
                print(line_out, flush=True)
                break
            check_docking = False
            stage_list += [[iteration, 'current']]
            write_stage(stage_list, stage_file)
            line_out = 'get_sample iteration: %d docking_done: %s' % (
                iteration, check_docking)
            print(line_out, flush=True)

            if auto_sub == 'shell':
                run_line = 'bash ' + sub_dock_script
                for i in range(num_sub):
                    subprocess.Popen(run_line.split())
                line_out = '%d sub_dock is generation' % num_sub
                print(line_out, flush=True)

            elif auto_sub == 'batch':
                run_line = 'sbatch ' + sub_dock_script
                for i in range(num_sub):
                    subprocess.Popen(run_line.split())
                line_out = '%d sub_dock_batch is submitted' % num_sub
                print(line_out, flush=True)

        if sleep_cycle is None:
            line_out = 'stop master: sleep_cycle is None '
            print(line_out, flush=True)
            break
        time.sleep(sleep_cycle)

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


def main():

    parser = argparse.ArgumentParser(description='Master of puppets')
    parser.register('action', 'extend', ExtendAction)
    parser.add_argument('--arg_file', type=open, required=False, default=None,
                        action=LoadFromConfig, help='argment file')
    parser.add_argument('--work_dir', type=str, required=False,
                        default='workflow', help='workflow directory')
    parser.add_argument('-s', '--smi_file_format', type=str, required=False,
                        default='pkl', help='pkl (default), txt, csv, tsv')
    parser.add_argument('--sleep_cycle', type=int, required=False, default=None,
                        help='sleep master for xx seconds per loop\n' +
                        'if None (default), escape without repeating the loop')
    parser.add_argument('--num_sample', type=int, required=False,
                        default='100000', help='number of samples')
    parser.add_argument('--file_size', type=int, required=False,
                        default='1000', help='smiles per file')
    parser.add_argument('--sampling_method', type=str, required=False,
                        default='sequential',
                        help='sequential, random, score, uncertainty')
    parser.add_argument('--mlsync', type=int, required=False, default=0,
                        help='0: sequential (alternative), 1: parallel')

    parser.add_argument('--auto_sub', type=str, required=False,
                        default=None,
                        help='shell (bash), batch (slurm) default is None')
    parser.add_argument('--sub_dock_script', type=str, required=False,
                        default='run_sub_dock.sh',
                        help='run_sub_dock.sh or slurm_sub_dock.sh')
    parser.add_argument('--num_sub', type=int, required=False, default=1,
                        help='number of docking subprocess')

    if len(sys.argv) < 2:
        parser.print_usage()
        sys.exit()
    args = parser.parse_args()

    work_dir = args.work_dir
    master_dir = work_dir + '/master'
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

    num_sample = args.num_sample
    sampling_method = args.sampling_method
    file_size = args.file_size
    sleep_cycle = args.sleep_cycle
    mlsync = args.mlsync

    auto_sub = args.auto_sub
    sub_dock_script = args.sub_dock_script
    num_sub = args.num_sub

    params_dict = dict()
    params_dict['work_dir'] = work_dir
    params_dict['master_dir'] = master_dir
    params_dict['todo_dir'] = todo_dir
    params_dict['current_dir'] = current_dir
    params_dict['done_dir'] = done_dir
    params_dict['field_separator'] = field_separator
    params_dict['smi_file_format'] = smi_file_format
    params_dict['num_sample'] = num_sample
    params_dict['sampling_method'] = sampling_method
    params_dict['file_size'] = file_size
    params_dict['sleep_cycle'] = sleep_cycle
    params_dict['mlsync'] = mlsync
    params_dict['auto_sub'] = auto_sub
    params_dict['sub_dock_script'] = sub_dock_script
    params_dict['num_sub'] = num_sub

    gapjil(params_dict)


if __name__ == "__main__":
    main()
