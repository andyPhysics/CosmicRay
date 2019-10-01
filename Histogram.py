import uproot
import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.pylab as plt
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Process Rooto Files')

parser.add_argument('-i',
                    dest = 'input_file',
                    help='This is the input root file')

parser.add_argument('-o',
                     dest = 'output_base',
                     help='Base of output files for histograms')

args = parser.parse_args()

path = '/data/ana/CosmicRay/IT73-IC79_3ySpectrumAndComposition/simulation/AllSim_merged'

file = path+'/'+args.input_file

f = uproot.open(file)

S125 = f['tinyTree']['s125'].array()

S125 = np.log10(S125)

Zenith= f['tinyTree']['zenith'].array()

EnergyLoss = zip(f['tinyTree']['eloss_1500'].array(),f['tinyTree']['eloss_1800'].array(),f['tinyTree']['eloss_2100'].array(),f['tinyTree']['eloss_2400'].array())

MeanEnergyLoss = [np.sum(i) for i in EnergyLoss]

MeanEnergyLoss = np.log10(MeanEnergyLoss)

HE_stoch_standard = f['tinyTree']['n_he_stoch'].array()

HE_stoch_strong = f['tinyTree']['n_he_stoch2'].array()

mass = f['tinyTree']['mass'].array()

Energy = f['tinyTree']['energy'].array()

Energy = np.log10(Energy)

fig, axes = plt.subplots(4, 2)

print(max(MeanEnergyLoss),min(MeanEnergyLoss))

axes[0,0].hist(S125,bins=10)
axes[0,0].set_title('Log10(S125)')
axes[0,1].hist(Zenith,bins=100)
axes[0,1].set_title('Zenith')
axes[1,0].hist(MeanEnergyLoss,bins=10)
axes[1,0].set_title('Log10(MeanEnergyLoss)')
axes[1,1].hist(HE_stoch_standard,bins=10)
axes[1,1].set_title('HE Stoch Standard')
axes[2,0].hist(HE_stoch_strong,bins=10)
axes[2,0].set_title('HE Stoch Standard')
axes[2,1].hist(mass,bins=100)
axes[2,1].set_title('Mass')
axes[3,0].hist(Energy,bins=100)
axes[3,0].set_title('Energy')
fig.tight_layout()
fig.savefig('Output'+args.output_base+'.png')

data = pd.DataFrame({'S125':S125,'Zenith':Zenith,'MeanEnergyLoss':MeanEnergyLoss,'HE_stoch_standard':HE_stoch_standard,'HE_stoch_strong':HE_stoch_strong,'mass':mass,'energy':Energy})

data.to_csv('output'+args.output_base+'.csv',index=True)


