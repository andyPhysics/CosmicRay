import pickle
import numpy as np
import uproot
from sklearn.metrics import confusion_matrix

def return_accuracy(predict1,predict2,validation):
    mass_predict1 = []
    for i in predict1:
        if i <2.5:
            mass_predict1.append(1)
        else:
            mass_predict1.append(4)
    cm=confusion_matrix(validation,mass_predict1)

    mass_predict4 = []
    for i in predict2:
        if i <2.5:
            mass_predict4.append(1)
        else:
            mass_predict4.append(4)
    cm1=confusion_matrix(validation,mass_predict4)
    
    return cm[0][0]/np.sum(cm[0]),cm1[1][1]/np.sum(cm1[1]),np.array(mass_predict1),np.array(mass_predict4)

file = uproot.open('/data/user/amedina/CosmicRay/Analysis/12632.root')
file2 = uproot.open('/data/user/amedina/CosmicRay/Analysis/12633.root')

cut_model = pickle.load((open('cut_model.sav','rb')))
proton_model =  pickle.load((open('proton_best_model.sav','rb')))
iron_model =  pickle.load((open('iron_best_model.sav','rb')))

S125_1 = np.log10(file['LaputopParams']['s125'].array())
S125_2 = np.log10(file2['LaputopParams']['s125'].array())
A1 = file['CurvatureOnlyParams']['A'].array()
A2 = file2['CurvatureOnlyParams']['A'].array()
D1 = file['CurvatureOnlyParams']['D'].array()
D2 = file2['CurvatureOnlyParams']['D'].array()
beta1 = file['LaputopParams']['beta'].array()
beta2 = file2['LaputopParams']['beta'].array()
zenith1 = file['Laputop']['zenith'].array()
zenith2 = file2['Laputop']['zenith'].array()
chi2_1 = file['CurvatureOnlyParams']['chi2_time'].array()
chi2_2 = file2['CurvatureOnlyParams']['chi2_time'].array()
energy1 = file['MCPrimary']['energy'].array()
energy2 = file2['MCPrimary']['energy'].array()
mass1 = [1 for i in range(len(energy1))]
mass2 = [4 for i in range(len(energy2))]

A = np.append(A1,A2)
D = np.append(D1,D2)
S125 = np.append(S125_1,S125_2)
beta = np.log10(np.append(beta1,beta2))
chi2 = np.log10(np.append(chi2_1,chi2_2))
zenith = np.cos(np.append(zenith1,zenith2))
zenith_new = np.append(zenith1,zenith2)*180/np.pi
mass = np.append(mass1,mass2)
energy = np.append(energy1,energy2)

mask = cut_model.predict(list(zip(S125,A,zenith,chi2)))

mask_new = []
for i in mask:
    if i < 0.5:
        mask_new.append(False)
    else:
        mask_new.append(True)

input_variable = np.array([np.append(i,j) for i,j in zip(D[mask_new],beta[mask_new])])
output = mass[mask_new]
energy = zenith_new[mask_new]


prediction1 = proton_model.predict(input_variable)
prediction2 = iron_model.predict(input_variable)
print(return_accuracy(prediction1,prediction2,output))

