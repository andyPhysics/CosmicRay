import numpy as np

def rotation(X,az,zen):
    R1 = np.array([[np.cos(az),-np.sin(az),0],
                [np.sin(az),np.cos(az),0],
                [0,0,1]])
    y_new = np.dot(R1,[0,1,0])
    R2_1 = np.cos(zen)*np.array([[1,0,0],[0,1,0],[0,0,1]])
    R2_2 = np.sin(zen)*(np.outer(np.cross(y_new,[1,0,0]),[1,0,0]) + np.outer(np.cross(y_new,[0,1,0]),[0,1,0]) + np.outer(np.cross(y_new,[0,0,1]),[0,0,1]))
    R2_3 = (1 - np.cos(zen)) * np.outer(y_new,y_new)
    R2 = R2_1+R2_2+R2_3
    X_prime  = np.dot(R2,np.dot(R1,X))
    return X_prime

def new_basis(az,zen):
    x_prime = rotation([1,0,0],az,zen)
    y_prime = rotation([0,1,0],az,zen)
    z_prime = rotation([0,0,1],az,zen)
    return x_prime,y_prime,z_prime

def new_vector(X,az,zen):
    x_prime,y_prime,z_prime = new_basis(az,zen)
    vector_x_prime = np.dot(x_prime,X)
    vector_y_prime = np.dot(y_prime,X)
    vector_z_prime = np.dot(z_prime,X)
    rho = ((vector_x_prime**2.0)+(vector_y_prime**2.0))**0.5
    return np.array([rho,vector_z_prime])
