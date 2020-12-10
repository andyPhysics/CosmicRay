import numpy as np
from sklearn.utils import resample
import matplotlib.pylab as plt
import matplotlib
from scipy.optimize import curve_fit

def weight_function(N,a,b):
    y = a * N + b
    return y

def new_df(df,mass):
    df['mass'] = [mass for i in range(len(df))]
    df['log_energy'] = np.log10(df['energy'])
    df['cos_zenith'] = np.cos(df['zenith'])
    df['log_energy_loss'] = np.log10(df['energy_loss'])
    df['log_m_r'] = np.log10(df['m_r'])
    df['log_s_r'] = np.log10(df['s_r'])
    df['log_charge'] = np.log10(df['charge'])
    df['log_A'] = np.log10(df['A'])
    m_o = df['m_o'].values
    m_125 = df['m_125'].values
    new = [(i-j) for i,j in zip(df['firstint'],df['Xo'])]
    df['new'] = new
    df['new_s125'] = [i + 125*j for i,j in zip(df['s_o'],df['s_r'])]
    df['delta_L'] = df['Xmax'] - df['firstint']
    energy_cut = []
    for i in df['log_energy']:
        if i > 7.0:
            energy_cut.append(True)
            continue
        elif np.random.random() < 0.5:
            energy_cut.append(True)
            continue
        energy_cut.append(False)
        
    fit = curve_fit(weight_function,df['N'][(df['waveform_weight']!=0)&(df['N']<75)],np.log10(df['waveform_weight'][(df['waveform_weight']!=0)&(df['N']<75)]))
    
    new_waveform_weight = []
    for i in range(len(df['waveform_weight'])):
        if (df['waveform_weight'].values[i] == 0)&(df['N'].values[i]<100):
            new_waveform_weight.append(10**weight_function(df['N'].values[i],fit[0][0],fit[0][1]))
            continue
        elif (df['waveform_weight'].values[i] == 0)&(df['N'].values[i]>=100):
            new_waveform_weight.append(1)
            continue
        new_waveform_weight.append(df['waveform_weight'].values[i])
    df['new_waveform_weight'] = new_waveform_weight
    
    df = df.loc[(df['Xmax']<900)&
                (df['Xmax']>400)&
                (np.array(energy_cut))&
                (df['ghRedChiSqr']<400)&
                (abs(df['max_check'])<20)&
                (df['S125']>0)&
                (df['fit_status_m']!=0)&
                (df['fit_status_m']!=2)&
                (df['m_chi2']>1e-4)&
                (df['log_m_r']>-5)&
                (df['log_m_r']<-1)&
                (df['m_125']<6)&
                (df['A']>1e-4)&
               (df['zenith']<45*(np.pi/180))]
    return df

def new_df_data(df):
    df['cos_zenith'] = np.cos(df['zenith'])
    df['log_energy_loss'] = np.log10(df['energy_loss'])
    df['log_m_r'] = np.log10(df['m_r'])
    df['log_s_r'] = np.log10(df['s_r'])
    df['log_charge'] = np.log10(df['charge'])
    df['log_A'] = np.log10(df['A'])
    m_o = df['m_o'].values
    m_125 = df['m_125'].values
    df['new_s125'] = [i + 125*j for i,j in zip(df['s_o'],df['s_r'])]
    df = df.loc[
                (df['S125']>0)&
                (df['fit_status_m']!=0)&
                (df['fit_status_m']!=2)&
                (df['m_chi2']>1e-4)&
                (df['log_m_r']>-5)&
                (df['log_m_r']<-1)&
                (df['m_125']<6)&
                (df['A']>1e-4)&
               (df['zenith']<45*(np.pi/180))]
    return df

def binning(x,y,bins,bins_input=[]):
    check1 = list(zip(x,y))
    if len(bins_input)!=0:
        hist1 = np.histogram(x,bins=bins_input)
    else:
        hist1 = np.histogram(x,bins=bins)
    mean = [[] for i in range(len(hist1[1][0:len(hist1[1])-1]))]
    mean_count = np.zeros(len(mean))
    index = range(len(mean))
    for i in check1:
        for j in index:
            if (i[0] > hist1[1][j] and i[0] < hist1[1][j+1]):
                mean[j].append( i[1])
                mean_count[j]+=1

    new_dist_10 = []
    new_dist_90 =  []
    check = []
    for j in mean:
        if len(j) == 0:
            check.append(False)
            continue
        else:
            check.append(True)
        dist10 = []
        dist90 = []
        for k in range(100):
            new_sample = resample(j)
            dist10.append(np.quantile(new_sample,0.25))
            dist90.append(np.quantile(new_sample,0.75))
        new_dist_10.append(dist10)
        new_dist_90.append(dist90)
    mean_overall = np.array([np.mean(i) for i in mean])[check]
    std_mean = np.array([np.std(i) for i in mean])[check]
    median = np.array([np.median(i) for i in mean])[check]
    ten = [np.mean(i) for i in new_dist_10]
    ninety = [np.mean(i) for i in new_dist_90]
    
    bins = np.array(hist1[1][0:len(hist1[1])-1])[check]

    return mean_overall,std_mean,bins,median,ten,ninety

def binning2(x,y,bins,confidence_interval = 1,bins_input=[]):
    check1 = list(zip(x,y))
    if len(bins_input)!=0:
        hist1 = np.histogram(x,bins=bins_input)
    else:
        hist1 = np.histogram(x,bins=bins)
    mean = [[] for i in range(len(hist1[1][0:len(hist1[1])-1]))]
    mean_count = np.zeros(len(mean))
    index = range(len(mean))
    for i in check1:
        for j in index:
            if (i[0] > hist1[1][j] and i[0] < hist1[1][j+1]):
                mean[j].append( i[1])
                mean_count[j]+=1

    new_dist_10 = []
    new_dist_90 =  []
    check = []
    for j in mean:
        if len(j) == 0:
            check.append(False)
            continue
        else:
            check.append(True)
        dist10 = []
        dist90 = []
        for k in range(100):
            new_sample = resample(j)
            if confidence_interval == 1:
                dist10.append(np.quantile(new_sample,0.25))
                dist90.append(np.quantile(new_sample,0.75))
            elif confidence_interval == 2:
                dist10.append(np.quantile(new_sample,0.10))
                dist90.append(np.quantile(new_sample,0.90))
        new_dist_10.append(dist10)
        new_dist_90.append(dist90)
    mean_overall = np.array([np.mean(i) for i in mean])[check]
    std_mean = np.array([np.std(i) for i in mean])[check]
    median = np.array([np.median(i) for i in mean])[check]
    ten = [np.mean(i) for i in new_dist_10]
    ninety = [np.mean(i) for i in new_dist_90]
    FWHM = np.array([2*(2*np.log(2))**0.5 * np.std(i) for i in mean])[check]
    
    bins = np.array(hist1[1][0:len(hist1[1])-1])[check]

    return mean,mean_overall,std_mean,bins,median,ten,ninety,FWHM

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
    
def cut_values(energy,mass_value = None,mass=[]):
    cut_energy = np.linspace(6.5,8,11)
    count = []
    for i in range(len(cut_energy[0:-1])):
        if (len(mass)==0)&(mass_value == None):
            count.append(np.sum((energy>cut_energy[i])&(energy<cut_energy[i+1])))
        else:
            count.append(np.sum((energy>cut_energy[i])&(energy<cut_energy[i+1])&(mass==mass_value)))
    return count

