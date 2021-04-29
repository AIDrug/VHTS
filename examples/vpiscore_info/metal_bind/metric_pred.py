import pandas as pd
import numpy as np

def cal_metric(df, num_active,num_data, value='Docking1'):
    df2 = df.sort_values(value, ascending=False)
    #print(df2)
    for num_select in range(100, df2.shape[0], 100):
        #num_select = 1000
        df3 = df2[0:num_select]
        ratio_a = (num_active/num_data)
        TP = (df3['activity']==True).sum()
        precision = TP/num_select
        print(precision/ratio_a, TP, num_select, num_active, num_data, precision, ratio_a)

file_name = 'workflow/master/docking.pkl'
#file_name = 'docking.pkl'

df1 = pd.read_pickle(file_name)
file_name2 = 'prediction.csv'
df2 = pd.read_csv(file_name2)
df = pd.merge(df1, df2[['MOL_ID', 'prediction']], on='MOL_ID')

num_data = df.shape[0]
print(num_data)
#num_active = 714
num_active = (df['activity']==True).sum()
print(num_active)
#cal_metric(df, num_active, value='Docking1')
#docking = df['Docking']
#pinfo = df['Pinfo']
#print(df[['MOL_ID', 'Pinfo']])
#pinfo2 = list()
#for i in range(num_data) :
#    pin = pinfo.loc[i]
#    pin2 = pin.clip(0,2)
#    a = np.argsort(-pin2)
#    pinfo2 += [pin2]

df_b1 = df[df['prediction']==True]
num_select = df_b1.shape[0]
print(num_select)
cal_metric(df_b1, num_active, num_data, value='VPIscore1')

#print(df_b1.loc[0])
df_b2 = df_b1.sort_values('VPIscore1', ascending=False)
df_b3 = df_b2[['MOL_IDX', 'MOL_ID', 'SMILES', 'activity', 'VPIscore_info1', 'VPIscore1', 'Docking1', 'Pinfo']]
file_out = 'docking_metal_rule.csv'
df_b3.to_csv(file_out, index=False)


