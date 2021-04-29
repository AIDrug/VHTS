import os
import pandas as pd
import numpy as np
from rdkit import Chem
from multiprocessing import Manager
from multiprocessing import Process
from multiprocessing import Queue
from openbabel import pybel
import subprocess


def ligand_preperation(smi, pH=7.4):
    run_line = 'obabel -:%s -osmi --neutralize' % (smi)
    result = subprocess.check_output(run_line.split(),
                                     stderr=subprocess.STDOUT,
                                     universal_newlines=True)
    for line in result.split('\n'):
        if line.startswith('1 molecule converted'):
            continue
        if line.startswith('0 molecule converted'):
            smi0 = smi
            break
        if len(line.strip()) < 1:
            continue
        smi0 = line.strip()
    try:
        m = pybel.readstring("smi", smi0)
        m.OBMol.AddHydrogens(False, True, pH)
        smi_p = m.write("smi").strip()
    except Exception:
        smi_p = smi0
    return smi_p


def creator(q, data, num_sub_proc):
    for d in data:
        idx = d[0]
        q.put((idx, d[1]))
    for i in range(0, num_sub_proc):
        q.put('DONE')


def worker(q, return_dict):
    smarts_pattern_list = [
                            '[#7;H0&R]~*~*~[#7;H0&R]',
#                            '[#7;H0&R]~*~*~[#8R0]',
                            '[#7;H0&R]~*~*~[#8R0;D1,D2&H1]',
                            '[#8R0;D1,D2&H1]~*~*~[#8R0;D1,D2&H1]',
                            '[#8R0;D1,D2&H1]~*~*~*~[#8R0;D1,D2&H1]',

#                            '[#8R0]~*~*~[#8R0]',
#                            '[#8R0]~*~*~*~[#8R0]',

#                            '[#7;H0&R]~*~*~*~[#8R0]',
#                            '[#8R0]~*~*~*~*~[#8R0]'
                            ]
    smarts_pattern_ring_list = ['[*r]1[*r][*r][*r]1',  '[*r]1[*r][*r][*r][*r]1',
                                '[*r]1[*r][*r][*r][*r][*r]1', '[*r]1[*r][*r][*r][*r][*r][*r]1',
                                '[*r]1[*r][*r][*r][*r][*r][*r][*r]1']
    smarts_list = list()
    for smarts_pattern in smarts_pattern_list:
        smarts = pybel.Smarts(smarts_pattern)
        smarts_list += [smarts]
    smarts_ring_list = list()
    for smarts_pattern in smarts_pattern_ring_list:
        smarts = pybel.Smarts(smarts_pattern)
        smarts_ring_list += [smarts]

    pout = 1000
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
        if idx % pout == pout-1:
            print("processing: ", idx+1, flush=True)

        prediction = False
        smi_p = ligand_preperation(smi, pH=7.4)
        try:
            m = pybel.readstring('smi', smi)
        except:
            return_dict[idx] = prediction
            continue

        m.addh()
        Nfeatures = 0
        features_m_list = list()
        for smarts in smarts_list:
            features = smarts.findall(m)
            Nfeatures += len(features)
            features_m_list += [features]
        if Nfeatures < 1:
            return_dict[idx] = prediction
            continue
        ring_list = list()
        for smarts in smarts_ring_list:
            ring_l0 = smarts.findall(m)
            for ring in ring_l0:
                ring_list += [set(ring)]

        Nfeatures_new = 0
        for ii, features in enumerate(features_m_list):
            for feature in features:
                check_subset = False
                if ii == 0:
                    aa = set((feature[0], feature[2]))
                    aa2 = set((feature[1], feature[3]))
                elif ii == 1 or ii == 2:
                    aa = set((feature[1], feature[2]))
                elif ii == 3 : # or ii == 4:
                    aa = set((feature[1], feature[3]))
                for ratoms in ring_list:
                    if aa.issubset(ratoms):
                        check_subset = True
                        break
                    if ii==0:
                        if aa2.issubset(ratoms):
                            check_subset = True
                            break

                if not check_subset:
                    Nfeatures_new += 1
        if Nfeatures_new > 0:
            prediction = True

        return_dict[idx] = prediction


def main():

    file_name = 'smiles_list.csv'
#    file_name = 'active.csv'

    df = pd.read_csv(file_name)
#    df = df[0:2000]
    num_data = df.shape[0]

    num_sub_proc = 30
    smiles_list = df[['MOL_IDX', 'SMILES']].values.tolist()
    data = list(enumerate(smiles_list))
    num_data = len(data)
    num_sub_proc = min(num_sub_proc, num_data)

    q1 = Queue()
    manager = Manager()
    return_dict = manager.dict()
    proc_master = Process(target=creator,
                          args=(q1, data, num_sub_proc))
    proc_master.start()

    procs = []
    for sub_id in range(0, num_sub_proc):
        proc = Process(target=worker, args=(q1, return_dict))
        procs.append(proc)
        proc.start()

    q1.close()
    q1.join_thread()
    proc_master.join()
    for proc in procs:
        proc.join()
    keys = sorted(return_dict.keys())

    prediction_list = list()
    for key in keys:
        prediction = return_dict[key]
        prediction_list += [prediction]
    df['prediction'] = prediction_list


    out_file = 'prediction.csv'
#    out_file = 'prediction_active.csv'

    df.to_csv(out_file, index=False)

    prediction = np.array(prediction_list, dtype=np.bool)
    activity = np.array(df['activity'], dtype=np.bool)
    NA = activity.sum()
    NP = prediction.sum()
    NTP = (activity*prediction).sum()
    print(NTP, NP, NA)
    print(NTP/NP, NTP/NA)


if __name__ == "__main__":
    main()
