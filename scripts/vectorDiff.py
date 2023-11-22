#!/usr/bin/env python3
#code author: James Beatty
#institution: The Ohio State Univesrity

import numpy as np
import numpy.linalg as LA
import math
import sys
import traceback

def tbexit():
    traceback.print_stack()
    sys.exit()

def isvectorarray(vec,dim=3):
    "is this an array of dim-vectors (True) or just a dim-vector (False)? Exit with traceback if neither."
    s = np.shape(vec) #gthis returns an empty tuple if vec is not an np.array
    if len(s) == 2 :
        if s[1] == dim:
            return True
        else:
            print("Error: isvectorarray called with 2d array, but not of (dim)-vectors. Exiting.")
            tbexit()
    elif len(s)==1 :
        if s[0]==dim : #have to separate these just in case s is empty
            return False
    else:
        print("Error: Bad argument to isvectorarray. Exiting.")
        tbexit()

def make_normal_basis(ref):
    tolerance = 0.01  #mimimum value of sin(theta) from z-axis for axis rotation
    xhat = np.array([1.,0.,0])
    zhat = np.array([0.,0.,1.])
    xp = np.cross(zhat,ref) 
    xpn=LA.norm(xp)
    if (xpn/LA.norm(ref) < tolerance):
        xp = xhat
    else:
        xp = xp/xpn
    yp = np.cross(ref,xp)
    return [xp,yp]
        
def space_angle_diff(ref,vec):
    "calculate 2d space angle differences of unit vector(s) from a reference unit vector"
    basis=make_normal_basis(ref)
    diff=np.cross(vec,ref)
    xpd=np.dot(diff,basis[0]) #positive xpd is toward zenith
    ypd=np.dot(diff,basis[1]) #positive ypd is clockwise around zenith, viewed from above.
    return np.arcsin(np.transpose(np.array([xpd,ypd])))

def normalize(vec,dim=3):
    if isvectorarray(vec,dim):
        retvec = np.array(vec)
        #need a copy here or else the argument gets overwritten
        for i in range((np.shape(vec))[0]):
            retvec[i,]=vec[i,]/LA.norm(vec[i,])
        return retvec
    else:
        return vec/LA.norm(vec)

def norm(vec,dim=3):
    if isvectorarray(vec,dim):
        retvec = np.array(np.zeros(np.shape(vec)[0]))
        #need a copy here or else the argument gets overwritten
        for i in range((np.shape(vec))[0]):
            retvec[i,]=LA.norm(vec[i,])
        return retvec
    else:
        return LA.norm(vec)

def gaussian_vector(len=10,sigma=1.,dim=2):
    rng=np.random.default_rng()
    retvec=np.array(np.zeros((len,dim)))
    for i in range(len):
        for j in range(dim):
            retvec[i,j]=rng.standard_normal() * sigma
    return retvec
    
# This is a test
def vdtest():
    refvec = normalize(np.array([1,0,1]))
    testvec = normalize(np.array([1,0.1,1])); 
    print(refvec,testvec,space_angle_diff(refvec,testvec))
    print()
    testvec2 = np.array([[1,0,1],[1,0.1,1],[1,0.2,1],[0,1,0],[-1,0,-1],[-1,-.1,-1],[1,0,0],[-1,0,0],[0,0,1],[1,1,0]])
    print(testvec2)
    testvec2=normalize(testvec2)
    result=space_angle_diff(refvec,testvec2)
    print(result)
    print(norm(result,2))
    random=gaussian_vector()
    print(random)
    samples=100000
    bigrandom=norm(gaussian_vector(samples,1/math.sqrt(2),2),2)
    #print(bigrandom)
    print(np.dot(bigrandom,bigrandom)/samples)

