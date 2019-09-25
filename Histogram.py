import uproot
import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.pylab as plt
import argparse

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

fig, axes = plt.subplots(3, 2)

print(max(MeanEnergyLoss),min(MeanEnergyLoss))

axes[0,0].hist(S125,bins=3)
axes[0,0].set_title('S125')
axes[0,1].hist(Zenith,bins=100)
axes[0,1].set_title('Zenith')
axes[1,0].hist(MeanEnergyLoss,bins=3)
axes[1,0].set_title('MeanEnergyLoss')
axes[1,1].hist(HE_stoch_standard,bins=10)
axes[1,1].set_title('HE Stoch Standard')
axes[2,0].hist(HE_stoch_strong,bins=10)
axes[2,0].set_title('HE Stoch Standard')
fig.tight_layout()
fig.savefig('Output'+args.output_base+'.png')



