#!/usr/bin/env python
import sys
import os
import argparse
from filelock import FileLock


def set_job_from_list(job_idx_list, list_dir):
    list_file = list_dir + '/list.txt'
    freeze_lock = FileLock('{}.lock'.format(list_file))
    with freeze_lock.acquire(timeout=30):
        if os.path.exists(list_file):
            with open(list_file, 'r') as fp:
                lines = fp.readlines()
        else:
            lines = list()
        with open(list_file, 'w') as fp:
            for job_idx in job_idx_list:
                line = job_idx + '\n'
                fp.write(line)
            for line in lines:
                fp.write(line)
    return


def remove_job_from_list(restart_list, list_dir):
    job_idx_list = list()
    list_file = list_dir + '/list.txt'
    freeze_lock = FileLock('{}.lock'.format(list_file))
    with freeze_lock.acquire(timeout=30):
        if os.path.exists(list_file):
            with open(list_file, 'r') as fp:
                lines = fp.readlines()
        with open(list_file, 'w') as fp:
            for line in lines:
                job_idx = line.strip()
                if restart_list != 'all':
                    if job_idx not in restart_list:
                        fp.write(line)
                    else:
                        job_idx_list += [job_idx]
                else:
                    job_idx_list += [job_idx]
    return job_idx_list


def move_current_to_todo(params_dict):
    restart_list = params_dict['restart_list']
    todo_dir = params_dict['todo_dir']
    current_dir = params_dict['current_dir']
    smi_file_format = params_dict['smi_file_format']
    job_idx_list = remove_job_from_list(restart_list, current_dir)
    print('restart job_idx list', job_idx_list, flush=True)
    for job_idx in job_idx_list:
        job_todo_file = todo_dir + '/' + job_idx + '.' + smi_file_format
        job_current_file = current_dir + '/' + job_idx + '.' + smi_file_format
        os.replace(job_current_file, job_todo_file)

    set_job_from_list(job_idx_list, todo_dir)

    return


def main():

    parser = argparse.ArgumentParser(description='worker for docking')
    parser.add_argument('--work_dir', type=str, required=False,
                        default='workflow', help='workflow directory')
    parser.add_argument('-s', '--smi_file_format', type=str, required=False,
                        default='pkl', help='pkl (default), txt, csv, tsv')
    parser.add_argument('--restart_list', type=str, required=True,
                        nargs='+',
                        help='example: --restart_list 0_0 0_1 0_2 ...' +
                        'or  --restart_list all')

    args = parser.parse_args()

    work_dir = args.work_dir
    todo_dir = work_dir + '/todo'
    current_dir = work_dir + '/current'
    done_dir = work_dir + '/done'
    smi_file_format = args.smi_file_format
    restart_list = args.restart_list
    if 'all' in restart_list:
       restart_list = 'all'
    params_dict = dict()
    params_dict['work_dir'] = work_dir
    params_dict['todo_dir'] = todo_dir
    params_dict['current_dir'] = current_dir
    params_dict['done_dir'] = done_dir
    params_dict['smi_file_format'] = smi_file_format
    params_dict['restart_list'] = restart_list

    move_current_to_todo(params_dict)


if __name__ == "__main__":
    main()
