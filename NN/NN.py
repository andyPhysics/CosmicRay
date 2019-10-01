import uproot
import numpy as np
import argparse
import pandas as pd

train_files = ["Helium_train.root","Proton_train.root","Iron_train.root","Oxygen_train.root"]
test_files = ["Helium_test.root","Proton_test.root","Iron_test.root","Oxygen_test.root"]
verification_files = ["Helium_verify.root","Proton_verify.root","Iron_verify.root","Oxygen_verify.root"]

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
        labels += y
        features += x
    labels = np.array(labels)
    features = np.array(features)
    return labels,features

    
train_labels,train_features = get_data(train_files)
test_labels,test_features = get_data(test_files)
print(test_features.shape)
import keras
from keras.layers import Dense, Dropout, Flatten, Input
from keras.models import Model

input_layer = Input(shape=(5,))

model = Dense(7,activation='tanh')(input_layer)

model = Dense(4,activation='tanh')(model)

predictions = Dense(2,activation='tanh')(model)

model = Model(inputs=input_layer,outputs=predictions)

opt=keras.optimizers.SGD()

model.compile(optimizer=opt , loss = 'mse')

model.fit(train_labels,train_features,
          epochs=100,
          validation_data = (test_labels,test_features))


