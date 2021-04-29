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

df = pd.read_pickle(file_name)
num_data = df.shape[0]
print(num_data)
#num_active = 714
num_active = (df['activity']==True).sum()
print(num_active)
cal_metric(df, num_active, num_data, value='VPIscore1')
#df_b2 = df_b1.sort_values('VPIscore1', ascending=False)

#file_out = 'docking_bond.csv'
#df_b2.to_csv(file_out, index=False)


