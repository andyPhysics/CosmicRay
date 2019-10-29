import numpy as np
from extract_data import *
import time
import matplotlib.pylab as plt
import multiprocessing as mp
from scipy.stats import chisquare

file_name = '/home/andy/Iron_all.npy'
loaded_dict = np.load(file_name,allow_pickle=True,encoding='latin1').item()

def concat_values(file_name,dict_key):
    x = np.load(file_name,allow_pickle=True,encoding='latin1').item()
    value = x[dict_key]      
    value_all = []
    for i in value:
        for j in i:
            value_all.append(j)
    return value_all
            
depth = concat_values(file_name,'depth')
depth = [i/1000.0 for i in depth]
run = concat_values(file_name,'run')
event = concat_values(file_name,'event')
E_plus_all = concat_values(file_name,'num_EPlus')
E_minus_all = concat_values(file_name,'num_EMinus')
E_all = [i+j for i,j in zip(E_plus_all,E_minus_all)]
E_all = [i/max(i) for i in E_all]
X_max = []
X_o = []
lambda_value = []
run_new = []
event_new = []
chi2_xmax = []

values = zip(run,event,depth,E_all)
count = 0
for i in values:
    print(count)
    
    try:
        output,depth_new,E_all_new = get_Xmax(i[2],i[3])
    except:
        continue
    
    chi2_xmax.append(chisquare(E_all_new,Gaisser_exp(depth_new,output[0],output[1],output[2],output[3]))[0])
    run_new.append(run[count])
    event_new.append(event[count])
    X_max.append(output[0]/output[1] + output[3])
    X_o.append(output[3])
    lambda_value.append(1/output[1])
    count+=1
    
new_dict = dict(run = run_new,
                event = event_new,
                X_max = X_max,
                X_o = X_o,
                lambda_value = lambda_value,
                chi2_xmax = chi2_xmax)

loaded_dict.update({'Gaisser_values':new_dict})
np.save('Iron_updated.npy',loaded_dict)



