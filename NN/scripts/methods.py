import numpy as np
from sklearn.utils import resample
import matplotlib.pylab as plt
import matplotlib
from scipy.optimize import curve_fit

def weight_function(N,a,b):
    y = a * N + b
    return y

def new_df(df, mass):
    add_basic_features(df, mass)
    add_derived_features(df)
    df = apply_energy_cut(df)
    fit = fit_waveform_weight(df)
    update_waveform_weight(df, fit)
    df = apply_additional_filters(df)
    
    return df

def new_df_data(df):
    add_data_features(df)
    df = apply_data_filters(df)
    
    return df

def add_basic_features(df, mass):
    df['mass'] = mass
    df['log_energy'] = np.log10(df['energy'])
    df['cos_zenith'] = np.cos(df['zenith'])
    df['log_muon_total'] = np.log10(df['muon_total'])
    df['log_energy_loss'] = np.log10(df['energy_loss'])
    df['log_m_r'] = np.log10(df['m_r'])
    df['log_s_r'] = np.log10(df['s_r'])
    df['log_charge'] = np.log10(df['charge'])
    df['log_A'] = np.log10(df['A'])

def add_derived_features(df):
    df['new'] = df['firstint'] - df['Xo']
    df['new_s125'] = df['s_o'] + 125 * df['s_r']
    df['delta_L'] = df['Xmax'] - df['firstint']

def apply_energy_cut(df):
    energy_cut_condition = (df['log_energy'] > 7.0) | (np.random.random(len(df)) < 0.5)
    return df[energy_cut_condition]

def fit_waveform_weight(df):
    fit = curve_fit(weight_function, df['N'][(df['waveform_weight'] != 0) & (df['N'] < 75)],
                    np.log10(df['waveform_weight'][(df['waveform_weight'] != 0) & (df['N'] < 75)]))
    return fit

def update_waveform_weight(df, fit):
    new_waveform_weight_condition = (df['waveform_weight'] == 0) & (df['N'] < 100)
    df.loc[new_waveform_weight_condition, 'new_waveform_weight'] = 10 ** weight_function(
        df.loc[new_waveform_weight_condition, 'N'], fit[0][0], fit[0][1])
    
    df.loc[(df['waveform_weight'] == 0) & (df['N'] >= 100), 'new_waveform_weight'] = 1

def apply_additional_filters(df):
    return df.loc[(df['Xmax'] < 900) & (df['Xmax'] > 400) &
                (df['ghRedChiSqr'] < 400) & (abs(df['max_check']) < 20) &
                (df['S125'] > 0) & (df['fit_status_m'].isin([0, 2])) &
                (df['m_chi2'] > 1e-4) & (-5 < df['log_m_r'] < -1) &
                (df['m_125'] < 6) & (df['A'] > 1e-4) &
                (df['zenith'] < 45 * (np.pi / 180))]

def add_data_features(df):
    df['cos_zenith'] = np.cos(df['zenith'])
    df['log_energy_loss'] = np.log10(df['energy_loss'])
    df['log_m_r'] = np.log10(df['m_r'])
    df['log_s_r'] = np.log10(df['s_r'])
    df['log_charge'] = np.log10(df['charge'])
    df['log_A'] = np.log10(df['A'])
    df['new_s125'] = df['s_o'] + 125 * df['s_r']

def apply_data_filters(df):
    return df.loc[
        (df['S125'] > 0) & (~df['fit_status_m'].isin([0, 2])) &
        (df['m_chi2'] > 1e-4) & (-5 < df['log_m_r'] < -1) &
        (df['m_125'] < 6) & (df['A'] > 1e-4) &
        (df['zenith'] < 45 * (np.pi / 180))]

def binning(x, y, bins, bins_input=None):
    bins_input = bins_input or bins
    hist1, bin_edges = np.histogram(x, bins=bins_input)
    indices = np.digitize(x, bin_edges[:-1])
    
    mean = [y[indices == i].tolist() for i in range(1, len(bin_edges))]
    
    new_dist_10 = [np.quantile(resample(j), 0.25) for j in mean]
    new_dist_90 = [np.quantile(resample(j), 0.75) for j in mean]
    
    mean_overall = np.array([np.mean(i) for i in mean if i]).tolist()
    std_mean = np.array([np.std(i) for i in mean if i]).tolist()
    median = np.array([np.median(i) for i in mean if i]).tolist()
    ten = np.array([np.mean(i) for i in new_dist_10])
    ninety = np.array([np.mean(i) for i in new_dist_90])
    
    return mean_overall, std_mean, bin_edges[:-1], median, ten, ninety

def binning2(x, y, bins, confidence_interval=1, bins_input=None):
    bins_input = bins_input or bins
    hist1, bin_edges = np.histogram(x, bins=bins_input)
    indices = np.digitize(x, bin_edges[:-1])
    
    mean = [y[indices == i].tolist() for i in range(1, len(bin_edges))]
    
    new_dist_10 = [np.quantile(resample(j), 0.10) for j in mean]
    new_dist_90 = [np.quantile(resample(j), 0.90) for j in mean]
    
    mean_overall = np.array([np.mean(i) for i in mean if i]).tolist()
    std_mean = np.array([np.std(i) for i in mean if i]).tolist()
    median = np.array([np.median(i) for i in mean if i]).tolist()
    ten = np.array([np.mean(i) for i in new_dist_10])
    ninety = np.array([np.mean(i) for i in new_dist_90])
    FWHM = np.array([2 * (2 * np.log(2))**0.5 * np.std(i) for i in mean if i]).tolist()
    
    return mean, mean_overall, std_mean, bin_edges[:-1], median, ten, ninety, FWHM

def plot_function(true,predicted,min_value,max_value,bins_hist,bins,name):
    bins_input = np.linspace(min_value,max_value,bins+1)
    percent_error = [100*(i-j)/j for i,j in zip(predicted,true)]
    mean_list,mean_overall,std_mean,bins,median,ten,ninety,FWHM = binning2(true,percent_error,bins=bins,confidence_interval = 2,bins_input=bins_input)
    mean_list1,mean_overall1,std_mean1,bins1,median1,ten1,ninety1,FWHM1 = binning2(true,predicted,bins=bins,confidence_interval = 2,bins_input=bins_input)
    
    fig, axs = plt.subplots(3, 2,figsize=(12,18))
    axs[0, 0].hist2d(true,predicted,bins=100,norm=matplotlib.colors.LogNorm())
    axs[0,0].plot(true,true,color='r',linestyle='--')
    axs[0, 0].set_title('Predicted vs True')
    axs[0,0].grid(True,linestyle='--')
    axs[0,1].scatter(bins1,FWHM1)
    axs[0, 1].set_title('FWHM vs True')
    axs[0,1].grid(True,linestyle='--')
    
    axs[1,0].hist2d(true,percent_error,bins=100,norm=matplotlib.colors.LogNorm())
    axs[1, 0].set_title('Percent Error vs True')
    axs[1,0].grid(True,linestyle='--')
    
    axs[1,1].scatter(bins,mean_overall)
    axs[1,1].set_title('Mean Error vs True')
    axs[1,1].grid(True,linestyle='--')
    
    axs[2,0].scatter(bins,ten,label='ten')
    axs[2,0].scatter(bins,ninety,label='ninety')
    axs[2,0].set_title('Confidence Intervals vs True')
    axs[2,0].legend()
    axs[2,0].grid(True,linestyle='--')
    plt.savefig(name+'1.png')
    
    count = 0
    fig2, axs2 = plt.subplots(4,4,figsize=(4*4,4*4))
    hist = np.histogram(percent_error,bins=bins_hist)
    for i in range(4):
        for j in range(4):
            if count+1>len(mean_list):
                break
            hist1 = axs2[i][j].hist(mean_list[count],bins=hist[1],edgecolor='white')
            axs2[i][j].set_title('%s - %s'%(bins_input[count],bins_input[count+1]))
            if len(mean_list[count]) != 0:
                axs2[i][j].vlines(ten[count],0,max(hist1[0]),color='red',linestyle='--')
                axs2[i][j].vlines(ninety[count],0,max(hist1[0]),color='red',linestyle='--')
                axs2[i][j].legend(['Mean: %.2f  \n Std: %.2f \n Ten: %.2f \n Ninety: %.2f'%(np.mean(mean_list[count]),
                                                                                            np.std(mean_list[count]),
                                                                                            ten[count],
                                                                                            ninety[count])],handlelength=0)
            axs2[i][j].grid(True,linestyle='--')
            count+=1
            
            
    plt.tight_layout()
    plt.savefig(name+'2.png')
    
def cut_values(energy, mass_value=None, mass=[]):
    cut_energy = np.linspace(6.5, 8, 11)
    count = []
    for i in range(len(cut_energy) - 1):
        mass_condition = (len(mass) == 0) & (mass_value is None) or (mass == mass_value)
        count.append(np.sum((energy > cut_energy[i]) & (energy < cut_energy[i + 1]) & mass_condition))
    return count

