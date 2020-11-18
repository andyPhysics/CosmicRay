#!/usr/bin/env python
# coding: utf-8


import tensorflow.keras as keras
from tensorflow.keras import initializers,regularizers
from tensorflow.keras.layers import Dense, Dropout, Flatten, Input, Concatenate, BatchNormalization, LeakyReLU
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.callbacks import EarlyStopping
import tensorflow.keras.backend as K
import math
import matplotlib.pylab as plt
import matplotlib
from random import random
from sklearn.utils import resample
import seaborn as sns
import corner
import sys,os
import pickle
from methods import * 
import matplotlib
import joblib


rs = int(sys.argv[1])


import numpy as np
import pandas as pd


from numpy.random import seed
seed(1)
from tensorflow.random import set_seed
set_seed(2)

#Load the data


iron = pd.read_csv('Iron.csv')
proton = pd.read_csv('Proton.csv')
helium = pd.read_csv('Helium.csv')
oxygen = pd.read_csv('Oxygen.csv')
proton2 = pd.read_csv('Proton2.csv')
iron2 =  pd.read_csv('Iron2.csv')
proton3 = pd.read_csv('Proton3.csv')
iron3 = pd.read_csv('Iron3.csv')
proton4 = pd.read_csv('Proton4.csv')
iron4 = pd.read_csv('Iron4.csv')

data = pd.read_csv('data.csv')


#Apply the cuts

iron = new_df(iron,4)
proton = new_df(proton,1)
helium = new_df(helium,2)
oxygen = new_df(oxygen,3)
proton2 = new_df(proton2,1)
iron2 = new_df(iron2,4)
proton3 = new_df(proton3,1)
iron3 = new_df(iron3,4)
proton4 = new_df(proton4,1)
iron4 = new_df(iron4,4)

#Create a master database

df = iron.append(proton)
df = df.append(helium)
df = df.append(oxygen)

#Ensuring that log_energy_loss is not nan

check = [math.isnan(i)==0 for i in df['log_energy_loss'].values]


df_coinc = df[check]

df_coinc = df_coinc.reset_index()

df_coinc.corr()[['Xmax','log_energy','mass','Xo']]


# ## Make train test splits and data

from sklearn.model_selection import train_test_split


test = df_coinc.sample(frac=0,random_state=rs)

output_variables = ['log_energy','Xmax','mass']
input_variables = ['cos_zenith','S125','log_energy_loss','he_stoch','he_stoch2','m_r','m_o','m_125','s_mean','s_std','A','waveform_weight']


test_y =  test[output_variables].values
test_X = test[input_variables].values

df_coinc = df_coinc.drop(test.index)
df_coinc = df_coinc.reset_index()


y = df_coinc[output_variables].values
X = df_coinc[input_variables].values

from sklearn.model_selection import KFold

kf = KFold(n_splits=5, random_state=42, shuffle=True)

count = 0
for train_index, test_index in kf.split(X):
    if count == rs:
        X_train, X_validation = X[train_index], X[test_index]
        y_train, y_validation = y[train_index], y[test_index]
    count += 1

#X_train, X_validation, y_train, y_validation = train_test_split(X, y, test_size=0.2, random_state=rs)


variable = list(zip(*X_train))
X_train1 = X_train[:,0:-1]
X_validation1 = X_validation[:,0:-1]


#Weight of models ensuring that masses are taken into consideration equallly

weight = []
for i in y_train:
    if i[0] == 1:
        weight.append(1)
    elif i[0] == 2:
        weight.append(len(proton)/len(helium))
    elif i[0] == 3:
        weight.append(len(proton)/len(oxygen))
    else:
        weight.append(len(proton)/len(iron))
weight=np.array(weight)


# ## Fit the Neural Network

model_best = 'model_coinc_best_%s.h5'%(rs)

best_model = keras.callbacks.ModelCheckpoint(model_best,
                                             monitor='val_loss',
                                             save_best_only=True,
                                             save_weights_only=False,
                                             mode='auto')

es = EarlyStopping(monitor='val_loss',patience=20)

input_layer = Input(shape=(len(X_train1[0]),))

model1 = BatchNormalization()(input_layer)

model1 = Dense(12,bias_regularizer=keras.regularizers.l1(1e-1))(model1)

model1 = LeakyReLU(alpha=0.2)(model1)

model2 = BatchNormalization()(model1)

model2 = Dense(7)(model2)

model2 = LeakyReLU(alpha=0.2)(model2)

model3 = Concatenate(axis=-1,activity_regularizer=keras.regularizers.l2(1e-5))([input_layer,model2])

model3 = Dropout(0.5)(model3)

model3 = BatchNormalization(renorm=True)(model3)

prediction1 = Dense(3,activation='linear',kernel_regularizer = keras.regularizers.l2(1e-5))(model3)

model = Model(inputs=input_layer,outputs=prediction1)

opt = keras.optimizers.Adam(decay=1e-5)

model.compile(optimizer=opt , loss = 'mse')


history = model.fit(X_train1,y_train[:,0:3],
                    epochs=200,
                    shuffle=True,
                    validation_data = (X_validation1,y_validation[:,0:3]),
                    callbacks=[best_model,es],
                    sample_weight=weight)


# ## Loss curve to observe how the network performs

best_model = load_model('model_coinc_best_%s.h5'%(rs))

# ## Fit the Decision Tree for Energy


from sklearn.ensemble import BaggingRegressor
from sklearn.tree import DecisionTreeRegressor


tree =BaggingRegressor(DecisionTreeRegressor(splitter='best',max_features='log2',random_state=42),n_estimators=400,bootstrap=True,random_state=42)


tree.fit(X_train1,y_train[:,0],sample_weight=weight)
joblib.dump(tree,'Energy_model_%s.pkl'%(rs))


predictions = best_model.predict(X_validation1)

predictions2 = tree.predict(X_validation1)


energy_predictions = predictions2
energy = np.array(list(zip(*y_validation))[0])
xmax = np.array(list(zip(*y_validation))[1])
xmax_predictions = np.array(list(zip(*predictions))[1])
mass = np.array(list(zip(*y_validation))[-1])
mass_predictions = np.array(list(zip(*predictions))[-1])


from scipy.optimize import curve_fit
def line_function(x,m,b):
    return m * x + b
def quadratic_function(x,m,b):
    y = b + m * x
    return y


value = [(i-j) for i,j in zip(xmax_predictions,xmax)]


#test_xmax = np.array(list(zip(*test_y)))[1]
#mass2 = np.array(list(zip(*test_y)))[2]


#check_predictions = best_model.predict(test_X[:,0:-1])
#check_predictions2 = tree.predict(test_X[:,0:-1])


from sklearn.linear_model import LinearRegression


line_model = LinearRegression()


line_model.fit(xmax_predictions.reshape(-1,1),value)
joblib.dump(line_model,'Xmax_bias_correction_%s.pkl'%(rs))

bias = line_model.predict(xmax_predictions.reshape(-1,1))
#bias2 = line_model.predict(test_xmax.reshape(-1,1))


xmax_predictions = np.array([(i-j) for i,j in zip(xmax_predictions,bias)])
#xmax_predictions2 = np.array([(i-j) for i,j in zip(list(zip(*check_predictions))[1],bias2)])
#mass_predictions2 = np.array(list(zip(*check_predictions))[-1])



check = (energy_predictions>=0)


predictions = predictions[check]


mass = mass[check]


energy = np.array(energy)



mass_check = np.array(list(zip(*y_validation))[-1])
#mass_check2 = np.array(list(zip(*test_y))[3])


#energy_true = list(zip(*test_y))[0]
energy_bins = np.linspace(6.5,8,11)


#xmax_true = np.array(list(zip(*test_y))[1])


mass_predictor =BaggingRegressor(DecisionTreeRegressor(splitter='best',random_state=42),random_state=42)


mass_predictor.fit(X_train1,y_train[:,-1],sample_weight=weight)
joblib.dump(mass_predictor,'mass_model_%s.pkl'%(rs))


mass_predictions = mass_predictor.predict(X_validation)


import scipy

def lnA_func(a,b,c,d):
    return 1 * a + 2 * b + 3 * c + 4 * d


from scipy.stats import chisquare


energy_bins = np.linspace(6.5,8,11)
count = 0

output_dict = {}

proton = []
iron = []
helium = []
oxygen = []
lnA = []

for m in range(5):
    for k in range(2):
        energy_check = []

        for j in energy_predictions:
            if (j>energy_bins[count])&(j<energy_bins[count+1]):
                energy_check.append(True)
            else:
                energy_check.append(False)

        energy_check = np.array(energy_check)
        
        
        lnA_list = []
        proton_list = []
        iron_list = []
        helium_list = []
        oxygen_list = []
        
        check_proton = (mass==1)&(energy_check)
        check_helium = (mass==2)&(energy_check)
        check_oxygen = (mass==3)&(energy_check)
        check_iron = (mass==4)&(energy_check)
        
        kde_proton = scipy.stats.gaussian_kde(xmax_predictions[check_proton])

        kde_helium= scipy.stats.gaussian_kde(xmax_predictions[check_helium])

        kde_oxygen = scipy.stats.gaussian_kde(xmax_predictions[check_oxygen])

        kde_iron = scipy.stats.gaussian_kde(xmax_predictions[check_iron])
            
        kde_proton_mass = scipy.stats.gaussian_kde(mass_predictions[check_proton])

        kde_helium_mass= scipy.stats.gaussian_kde(mass_predictions[check_helium])

        kde_oxygen_mass = scipy.stats.gaussian_kde(mass_predictions[check_oxygen])

        kde_iron_mass = scipy.stats.gaussian_kde(mass_predictions[check_iron])
            
        def my_function(x,a,b,c,d):
            y1 = kde_proton.evaluate(x)
            y2 = kde_helium.evaluate(x)
            y3 = kde_oxygen.evaluate(x)
            y4 = kde_iron.evaluate(x)
            y = a * y1 + b * y2 + c * y3 + d * y4
            return y
        
        def my_function2(x,a,b,c,d):
            y1 = kde_proton_mass.evaluate(x) 
            y2 = kde_helium_mass.evaluate(x) 
            y3 = kde_oxygen_mass.evaluate(x)
            y4 = kde_iron_mass.evaluate(x)
            y = a * y1 + b * y2 + c * y3 + d * y4
            return y
       
        error = np.array([(i-j) for i,j in zip(xmax_predictions,y_validation[:,1])])
    
        output_dict['%s - %s'%(energy_bins[count],energy_bins[count+1])] = {}
        output_dict['%s - %s'%(energy_bins[count],energy_bins[count+1])]['proton'] = xmax_predictions[check_proton]
        output_dict['%s - %s'%(energy_bins[count],energy_bins[count+1])]['helium'] = xmax_predictions[check_helium]
        output_dict['%s - %s'%(energy_bins[count],energy_bins[count+1])]['oxygen'] = xmax_predictions[check_oxygen]
        output_dict['%s - %s'%(energy_bins[count],energy_bins[count+1])]['iron'] = xmax_predictions[check_iron]
        
        output_dict['%s - %s'%(energy_bins[count],energy_bins[count+1])]['proton_error'] = error[check_proton]
        output_dict['%s - %s'%(energy_bins[count],energy_bins[count+1])]['helium_error'] = error[check_helium]
        output_dict['%s - %s'%(energy_bins[count],energy_bins[count+1])]['oxygen_error'] = error[check_oxygen]
        output_dict['%s - %s'%(energy_bins[count],energy_bins[count+1])]['iron_error'] = error[check_iron]
        
        output_dict['%s - %s'%(energy_bins[count],energy_bins[count+1])]['proton_mass'] = mass_predictions[check_proton]
        output_dict['%s - %s'%(energy_bins[count],energy_bins[count+1])]['helium_mass'] = mass_predictions[check_helium]
        output_dict['%s - %s'%(energy_bins[count],energy_bins[count+1])]['oxygen_mass'] = mass_predictions[check_oxygen]
        output_dict['%s - %s'%(energy_bins[count],energy_bins[count+1])]['iron_mass'] = mass_predictions[check_iron]
        
        count += 1
        

np.save('kde_plots_%s'%(rs),output_dict)

