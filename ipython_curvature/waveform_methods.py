import numpy as np
import matplotlib.pylab as plt
import plotly.graph_objects as go
from plotly.offline import plot
from ipywidgets import widgets
import uproot
from icecube.icetray import *


def process_file(run,event,input_file,geometry):
    #-----------upload the file-------------
    input_file = np.load(input_file,allow_pickle=True,encoding='latin1').item()
    geometry = np.load(geometry,allow_pickle=True,encoding='latin1').item()
    #------------define constants----------------
    c = .299 #m/ns
    #----------This is the output dictionary-------------------
    output = {}
    #------------This is the geometry--------------------------
    output['icetop_geometry'] = {}
    geometry_keys = geometry['geometry'].keys()
    string = []
    dom = []
    icetop = [61,62,63,64]
    
    for i in geometry_keys:
        if i.om not in icetop:
            continue
        string.append(i.string)
        dom.append(i.om)

    for i in string:
        for j in dom:
            output['icetop_geometry'][OMKey(i,j,0)] = {}
            output['icetop_geometry'][OMKey(i,j,0)]['position'] = []
            output['icetop_geometry'][OMKey(i,j,0)]['angles'] = []
                
    for i in output['icetop_geometry'].keys():
        output['icetop_geometry'][i]['position'].append(geometry['geometry'][i].x)
        output['icetop_geometry'][i]['position'].append(geometry['geometry'][i].y)
        output['icetop_geometry'][i]['position'].append(geometry['geometry'][i].z)
        output['icetop_geometry'][i]['angles'].append(geometry['geometry'][i].theta)
        output['icetop_geometry'][i]['angles'].append(geometry['geometry'][i].phi)
        
    
    output[run] = {}
    output[run][event] = {}
    
        
    #-----------------update output to include important variables of the different events-----------------
    xc = input_file[run][event]['x']
    yc = input_file[run][event]['y']
    zpos = input_file[run][event]['z']
    tc = input_file[run][event]['time']
    azimuth = input_file[run][event]['azimuth']
    zenith = input_file[run][event]['zenith']
    energy = input_file[run][event]['energy']

    output[run][event]['event_info'] = {}
    output[run][event]['event_info']['xc'] = xc
    output[run][event]['event_info']['yc'] = yc
    output[run][event]['event_info']['zpos'] = zpos
    output[run][event]['event_info']['tc'] = tc
    output[run][event]['event_info']['energy'] = energy
    output[run][event]['event_info']['azimuth'] = azimuth
    output[run][event]['event_info']['zenith'] = zenith
    
    
    output[run][event]['SLC_Info'] = {}
    output[run][event]['SLC_Info']['charge'] = {}
    output[run][event]['SLC_Info']['time'] = {}
    output[run][event]['SLC_Info']['r'] = {}
    for DOM in input_file[run][event]['SLC_Info']['charge'].keys():
        output[run][event]['SLC_Info']['charge'][DOM] = input_file[run][event]['SLC_Info']['charge'][DOM]
        output[run][event]['SLC_Info']['time'][DOM] = input_file[run][event]['SLC_Info']['time'][DOM]
        
        position = output['icetop_geometry'][DOM]['position']
        x_dom = position[0]
        y_dom = position[1]
        z_dom = position[2]
        
        output[run][event]['SLC_Info']['r'][DOM] = ((xc-x_dom)**2.0 + (yc-y_dom)**2.0)**0.5
        
         
                    
    for DOM in input_file[run][event]['waveforms'].keys():
        x = input_file[run][event]['waveforms'][DOM]
        start_time = input_file[run][event]['waveform_time'][DOM]
        z = input_file[run][event]['waveform_binwidth'][DOM]
        
        position = output['icetop_geometry'][DOM]['position']
        x_dom = position[0]
        y_dom = position[1]
        z_dom = position[2]
            
        phi_dom = output['icetop_geometry'][DOM]['angles'][0]
        theta_dom = output['icetop_geometry'][DOM]['angles'][1]
        
        time_list1 = [i*z + start_time - tc for i in range(len(x))]
        
        #-----------Calculate CDF------------
                        
        output[run][event][DOM] = {}            
        count1 = 0
            
        output[run][event][DOM]['waveform'] = x
        output[run][event][DOM]['time'] = time_list1
        output[run][event][DOM]['CDF'] = []
                
            
        for i in range(len(x)):
            if count1 == 0:
                output[run][event][DOM]['CDF'].append(x[i])
            else:
                output[run][event][DOM]['CDF'].append(x[i]+output[run][event][DOM]['CDF'][i-1])
            count1+=1
        #-------------------------------------  
        DOM_string = str(DOM)    
        output[run][event][DOM]['time_10'] = input_file[run][event]['waveform_info'][DOM_string]['Time_10']
        output[run][event][DOM]['time_50'] = input_file[run][event]['waveform_info'][DOM_string]['Time_50']
        output[run][event][DOM]['time_90'] = input_file[run][event]['waveform_info'][DOM_string]['Time_90']
        output[run][event][DOM]['charge'] = input_file[run][event]['waveform_info'][DOM_string]['Charge_PE']
        output[run][event][DOM]['r'] = ((xc-x_dom)**2.0 + (yc-y_dom)**2.0)**0.5
        
    for DOM in input_file[run][event]['All_pulses'].keys():
        output[run][event][DOM]['all_time'] = input_file[run][event]['All_pulses'][DOM]['time']
        output[run][event][DOM]['all_charge'] = input_file[run][event]['All_pulses'][DOM]['charge']
        output[run][event][DOM]['all_radius'] = input_file[run][event]['All_pulses'][DOM]['radius']
            
    return output

def enable_plotly_in_cell():
  import IPython
  from plotly.offline import init_notebook_mode
  display(IPython.core.display.HTML('''
        <script src="/static/components/requirejs/require.js"></script>
  '''))
  init_notebook_mode(connected=True)

def get_event(loaded_file,run,event,length_of_arrow,output_html,auto_open):

    enable_plotly_in_cell

    geometry_x = []
    geometry_y = []
    geometry_z = []
    relative_azimuth = []
    relative_zenith = []
    DOM_key = []
    color_scale = []
    for dom in loaded_file['icetop_geometry'].keys():
 
        DOM_key.append('%s'%(dom))
        geometry_x.append(loaded_file['icetop_geometry'][dom]['position'][0])
        geometry_y.append(loaded_file['icetop_geometry'][dom]['position'][1])
        geometry_z.append(loaded_file['icetop_geometry'][dom]['position'][2])
        relative_azimuth.append('relative_azimuth: %.2f  %.2f'%(loaded_file[run][event]['event_info']['azimuth'],loaded_file['icetop_geometry'][dom]['angles'][1]))
    
        if dom in loaded_file[run][event].keys():
            #color_scale.append(np.log10(loaded_file[run][event][dom]['time_90']-loaded_file[run][event][dom]['time_10']))
            color_scale.append(np.log10(loaded_file[run][event][dom]['charge']))
        else:
            color_scale.append('white')
        
    text_info = [i+'<br>'+j for i,j in zip(DOM_key,relative_azimuth)]

    length = length_of_arrow
    azimuth = loaded_file[run][event]['event_info']['azimuth']
    zenith = loaded_file[run][event]['event_info']['zenith']

    vector_x = loaded_file[run][event]['event_info']['xc']
    vector_x_1 = loaded_file[run][event]['event_info']['xc'] + length * np.sin(zenith) * np.cos(azimuth)

    vector_y = loaded_file[run][event]['event_info']['yc']
    vector_y_1 = loaded_file[run][event]['event_info']['yc'] + length * np.sin(zenith) * np.sin(azimuth)

    vector_z = loaded_file[run][event]['event_info']['zpos']
    vector_z_1 = loaded_file[run][event]['event_info']['zpos'] + length * np.cos(zenith)

    layout = go.Layout(
            scene = dict(
            xaxis = dict(nticks=10, range=[-700,700],),
                         yaxis = dict(nticks=10, range=[-700,700],),
                         zaxis = dict(nticks=10, range=[1940,3340],),),
            width=700,
            margin=dict(r=20, l=10, b=10, t=10))            
 
    fig = go.FigureWidget(layout=layout)            
            
    geometry = go.Scatter3d(x=geometry_x, y=geometry_y, z=geometry_z,
                            mode='markers',
                            marker=dict(size=3,opacity=0.5,
                            color=color_scale,colorscale='Inferno',
                            showscale=True),
                            name = 'DOMs',
                            text = text_info)

    fig.add_trace(geometry)

    vector = fig.add_trace(go.Scatter3d( x = [vector_x,vector_x_1],
                           y = [vector_y,vector_y_1],
                           z = [vector_z,vector_z_1],
                          mode = 'lines+markers',
                           marker = dict( size = 5,
                                          color = "rgb(84,48,5)"),
                           line = dict( color = "rgb(84,48,5)",
                                        width = 6),
                         name = "Event",
                         text = 'Energy: ' + str(np.log10(loaded_file[run][event]['event_info']['energy']))))

            
    plot(fig,filename=output_html,auto_open=auto_open)


def waveform_output(loaded_file,run,event,output_html,auto_open,CDF):
    enable_plotly_in_cell       
    gain = np.load('DOM_gain.npy',allow_pickle=True,encoding='latin1').item()
    keys = list(loaded_file[run][event].keys())[2:]
    CDF = CDF
    if CDF:
        output_value = 'CDF'
    else:
        output_value = 'waveform'

    plots = []
    plots2 = []
    plots3 = []
    plots4 = []
    update_menu = []
    update_menu2 = []

    count = 0
    check = None
    for i in keys:
        string = str(i).split(',') #This is the tank number
        test2 = string[0].split('(')[1] #This is the string number
        visible = [False]*len(keys)*2
        visible[count] = True
        visible[len(keys) + count] = True
        visible2 = [False]*len(keys)*2
        visible_out = visible + visible2
        
        plots.append(go.Scatter(x=list(loaded_file[run][event][i]['time']),
                                y=list(loaded_file[run][event][i][output_value]),
                               mode='markers',
                               name = str(i)+' Gain:%s'%(gain[i]),
                               visible=False))
        plots2.append(go.Scatter(x=list(loaded_file[run][event][i]['time']),
                                y=list(loaded_file[run][event][i][output_value]),
                                name = str(i),
                                 visible=False,
                                showlegend=False))
        update_menu.append(dict(
                        args=[{'visible': visible_out}],
                        label=str(i) + " CDF",
                        method="restyle"
                    ))
        count += 1
        
    count = 0
    check = None
    for i in keys:
        string = str(i).split(',') #This is the tank number
        test2 = string[0].split('(')[1] #This is the string number
        check_list = [61,62]
        check_list2 = [63,64]
        if check != test2:
            check = test2
            visible = [False]*len(keys)*2
            visible2 = [False]*len(keys)*2
            visible2[count] = True
            visible2[len(keys) + count] = True
        elif check == test2:
            visible2[count] = True
            visible2[len(keys) + count] = True
        
        
        visible_out = visible + visible2
        
        plots3.append(go.Scatter(x=list(loaded_file[run][event][i]['time']),
                                y=list(loaded_file[run][event][i]['waveform']),
                                mode='markers',
                                visible= False,
                                name = str(i)+' Gain:%s'%(gain[i])))
        plots4.append(go.Scatter(x=list(loaded_file[run][event][i]['time']),
                                y=list(loaded_file[run][event][i]['waveform']),
                                name = str(i),
                                visible=False,
                                showlegend=False))
        update_menu2.append(dict(
                        args=[{'visible': visible_out}],
                        label=str('DOM_'+test2) + ' waveform',
                        method="restyle"
                    ))
        count += 1
        
    plot_output = plots + plots2 + plots3 + plots4

    updatemenus = [go.layout.Updatemenu(buttons=update_menu,
                                            direction="down",
                                            pad={"r": 10, "t": 10},
                                            showactive=True,
                                            x=0.1,
                                            xanchor="left",
                                            y=1.1,
                                            yanchor="top"
                                            ),
                   go.layout.Updatemenu(buttons=update_menu2,
                                            direction="down",
                                            pad={"r": 10, "t": 10},
                                            showactive=True,
                                            x=0.1,
                                            xanchor="right",
                                            y=1.1,
                                            yanchor="top"
                                            )
                        ]  
    
        
    layout = go.Layout(
        xaxis_title = 'Time [ns]',
        yaxis_title = output_value + ' [mV]',
        updatemenus=updatemenus
    )


    fig = go.Figure(data=plot_output, layout=layout)
    plot(fig,filename=output_html,auto_open=auto_open)

def binning(x,y,bins):
    check1 = list(zip(x,y))
    hist1 = np.histogram(x,bins=bins)
    mean = [[] for i in range(len(hist1[1][0:len(hist1[1])-1]))]
    mean_count = np.zeros(len(mean))
    index = range(len(mean))
    for i in check1:
        for j in index:
            if (i[0] > hist1[1][j] and i[0] < hist1[1][j+1]):
                mean[j].append( i[1])
                mean_count[j]+=1

    mean_overall = [np.mean(i) for i in mean]
    std_mean = [np.std(i) for i in mean]
    bins = hist1[1][0:len(hist1[1])-1]
    return mean_overall,std_mean,bins
