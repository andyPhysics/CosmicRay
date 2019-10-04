import uproot
import numpy as np
import argparse
import pandas as pd
import keras.backend as K
from tensorflow.python.framework import ops
from tensorflow.python.ops import gen_math_ops as math_ops


verification_files = ["Helium_verify.root","Proton_verify.root","Iron_verify.root","Oxygen_verify.root"]

from sklearn.preprocessing import minmax_scale
from keras.models import load_model

#----------------------------------------------------------------------------
def get_data(input_file_list):
    path = '/data/ana/CosmicRay/IT73-IC79_3ySpectrumAndComposition/simulation/AllSim_merged'
    filelist = input_file_list
    labels = []
    features = []

    for i in filelist:
        f = uproot.open(path+'/'+i)
        Energy = f['tinyTree']['energy'].array()
        Energy = np.log10(Energy)
        Mass = f['tinyTree']['mass'].array()
        Mass = np.log(Mass)
        S125 = f['tinyTree']['s125'].array()
        S125 = np.log10(S125)
        Zenith= f['tinyTree']['zenith'].array()
        EnergyLoss = zip(f['tinyTree']['eloss_1500'].array(),f['tinyTree']['eloss_1800'].array(),f['tinyTree']['eloss_2100'].array(),f['tinyTree']['eloss_2400'].array())
        MeanEnergyLoss = [np.mean(i) for i in EnergyLoss]
        MeanEnergyLoss = np.log10(MeanEnergyLoss)
        HE_stoch_standard = f['tinyTree']['n_he_stoch'].array()
        HE_stoch_strong = f['tinyTree']['n_he_stoch2'].array()
       

        x = zip(Energy,Mass)
        y = zip(S125,np.cos(Zenith),MeanEnergyLoss,HE_stoch_standard,HE_stoch_strong)
        features += y
        labels += x
#    features = [minmax_scale(i,feature_range=(-1,1)) for i in features]
    features = np.array(features)
    labels = np.array(labels)
    return labels,features


def custom_loss(ytrue,ypred):
    y_pred1 = ops.convert_to_tensor(ypred[0])
    y_pred2 = ops.convert_to_tensor(ypred[1])

    y_true1 = math_ops.cast(ytrue[0], ypred[0].dtype)
    y_true2 = math_ops.cast(ytrue[1], ypred[1].dtype)

    return K.mean(math_ops.square(y_pred1 - y_true1), axis=-1)+10.0*K.mean(math_ops.square(y_pred2 - y_true2), axis=-1)


    
labels,features = get_data(verification_files)

model = load_model('NN_best.h5',custom_objects={'custom_loss':custom_loss})

label_1 = []

for i in verification_files:
    x,y = get_data([i])
    labels_pred = model.predict(y)
    file_output = zip(x,labels_pred)
    label_1.append(file_output)

output_labels = model.predict(features)

output = zip(labels,output_labels)

np.save('Output.npy',output)
np.save('Output_split.npy',label_1)
