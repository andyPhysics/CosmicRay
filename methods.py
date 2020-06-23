from vectorDiff import *
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from functools import partial

from I3Tray import *
from icecube import icetray, dataclasses, dataio, toprec, recclasses, frame_object_diff,simclasses,WaveCalibrator,tpx
load('millipede')
load('stochastics')
from icecube.icetray import I3Module, OMKey



c = .299

def get_delta(r,a,b,sigma):
    b_start = 19.51
    sigma_start = 83.5
    delta = a * r**2.0 #+ b_start * (1-np.exp(-(r**2.0)/(2*(sigma_start**2.0))))
    return delta

def get_n(x,y,z):
    x,y,z = np.array([x,y,z])/(x**2.0+y**2.0+z**2.0)**0.5
    return [x,y,z]

def get_n_cos(zen,az):
    x = -np.sin(zen) * np.cos(az)
    y = -np.sin(zen) * np.sin(az)
    z = -np.cos(zen)
    return x,y,z

def get_t(X,x1,y1,z1,a,b,sigma,xc,yc,zc,tc):
    x,y,z,r = X
    n = get_n(x1,y1,z1)
    x_new = np.array([(i-xc)*n[0] for i in x])
    y_new = np.array([(i-yc)*n[1] for i in y])
    z_new = np.array([(i-zc)*n[2] for i in z])
    new = x_new + y_new + z_new
    tc = np.array([tc for i in range(len(z_new))])
    t = tc + (1/c)*new + get_delta(r,a,b,sigma)
    return t

def get_ang_diff(x1,y1,z1,x2,y2,z2):
    n1 = np.array(get_n(x1,y1,z1))
    n2 = np.array(get_n(x2,y2,z2))
    if np.dot(n1,n2)>1:
        value = 1
        angular_diff = np.arccos(value)*180/np.pi
    else:
        angular_diff = np.arccos(np.dot(n1,n2))*180/np.pi
    return angular_diff

def magnitude_spherical(theta,d_theta,d_phi):
    dl = (d_theta)**2.0 + (np.sin(theta)**2.0)*(d_phi)**2.0
    return dl

def get_fit(frame,tol):    
    true_zen = frame['MCPrimary'].dir.zenith
    true_az = frame['MCPrimary'].dir.azimuth
    
    true_x,true_y,true_z = get_n_cos(true_zen,true_az)

    tc = frame['ShowerCOG'].time
    xc = frame['ShowerCOG'].pos.x
    yc = frame['ShowerCOG'].pos.y
    zc = frame['ShowerCOG'].pos.z
    zen_shower = frame['ShowerPlane'].dir.zenith
    az_shower = frame['ShowerPlane'].dir.azimuth
    
    x_start,y_start,z_start = get_n_cos(zen_shower,az_shower)
    a_start = 4.823*10**-4
    b_start = 19.51
    sigma_start = 83.5
    
    x_o = []
    y_o = []
    z_o = []
    r_o = []
    t0_o = []
    for omkey in frame['WaveformInfo'].keys():
        dom = int(omkey.split(',')[1])
        string = int(omkey.split(',')[0].split('(')[1])
        position = frame['I3Geometry'].omgeo[OMKey(string,dom)].position
        x_o.append(position.x)
        y_o.append(position.y)
        z_o.append(position.z)
        r_o.append(((position.x)**2.0 + (position.y)**2.0 + (position.z)**2.0)**0.5)
        t0_o.append(frame['WaveformInfo'][omkey]['t_0'])


    x_o = np.array(x_o)
    y_o = np.array(y_o)
    z_o = np.array(z_o)
    r_o = np.array(r_o)
    t0_o = np.array(t0_o)
    
    get_t_new = partial(get_t,tc=tc,xc=xc,yc=yc,zc=zc)
    
    check = [True for i in range(len(event1))]
    remove_index_list = []
    
    fit_original = curve_fit(get_t_new,(x_o,y_o,z_o,r_o),t0_o,
                    p0=[x_start,y_start,z_start,a_start,b_start,sigma_start],
                    bounds=((-1,-1,-1,0,0,0),(1,1,0,2e-3,200,500)),
                    maxfev=3000,
                    absolute_sigma=True)

    vector_list = []
    fit_list = []
    mag_list = [0]
    x,y,z = get_n(fit_original[0][0],fit_original[0][1],fit_original[0][2])
    vector_list.append(np.array([x,y,z]))
    fit_list.append(fit_original)

    pass_value = False

    count = 0
    
    while True:
        if count == 0:
            check1 = np.copy(check)
        else:
            check1[remove_index] = False

        fit_check = []
        removed = []
        vector_check_list = []
        
        for i in np.array(range(len(check1))):
            check2 = np.copy(check1)
            check2[i] = False

            if (i in remove_index_list):
                continue

            x = x_o[check2]
            y = y_o[check2]
            z = z_o[check2]
            r = r_o[check2]
            t0 = t0_o[check2]
    
            fit = curve_fit(get_t_new,(x,y,z,r),t0,
                            p0=[fit_list[-1][0][0],fit_list[-1][0][1],fit_list[-1][0][2],fit_list[-1][0][3],fit_list[-1][0][4],fit_list[-1][0][5]],
                            bounds=((-1,-1,-1,0,0,0),(1,1,0,2e-3,200,500)),
                            maxfev=5000,
                            absolute_sigma=True)
            fit_check.append(fit)
            x,y,z = get_n(fit[0][0],fit[0][1],fit[0][2])
            vector_check_list.append([x,y,z])
            removed.append(i)
        
        x_o,y_o,z_o = get_n(fit_list[-1][0][0],fit_list[-1][0][1],fit_list[-1][0][2])
        
        theta = np.arctan(((x_o**2.0+y_o**2.0)**0.5)/z_o)
        
        output = space_angle_diff(vector_list[-1],vector_check_list)
        
        mag = [magnitude_spherical(theta,d_theta,d_phi)**0.5 for d_theta,d_phi in output]
        
        remove_index = removed[np.argmax(mag)]
        remove_index_list.append(remove_index)
        
        fit_list.append(fit_check[np.argmax(mag)])
        
        mag_list.append(mag[np.argmax(mag)])
        
        vector1 = np.sum(vector_check_list,axis=0)/len(vector_check_list)
        vector2 = np.array([x_o,y_o,z_o])
        bias = (len(vector_check_list)-1)*(vector1-vector2)
        
        if abs(mag_list[-1]-mag_list[-2])<tol:
            check1[remove_index] = False
            final_fit = fit_check[np.argmax(mag)]
            x_final,y_final,z_final = get_n(final_fit[0][0],final_fit[0][1],final_fit[0][2])
            ang_diff_final = get_ang_diff(x_final,y_final,z_final,true_x,true_y,true_z)
            break
            
        
        
    return final_fit, check1

class New_fit(I3Module):
    def __init__(self, context):
        I3Module.__init__(self, context)

    def Physics(self, frame):
        fit, mask = get_fit(fit,0.5)
        mask1 = dataclasses.I3RecoPulseSeriesMapMask(self.frame, 'LaputopHLCVEM')
        count = 0
        for omkey in frame['WaveformInfo'].keys():
            dom = int(omkey.split(',')[1])
            string = int(omkey.split(',')[0].split('(')[1])
            key = OMKey(string,dom)
            mask1.set(key, 1, mask[count])


        self.PushFrame(frame)
