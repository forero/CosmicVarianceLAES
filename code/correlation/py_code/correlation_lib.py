import numpy as np


def DD_histogram(X,Y,distance,th_min,th_max,theta_bins=10,cat_number=1):

    n_points=len(X)
    d_arr=[]
    for i in range(n_points):
        for j in range(i,n_points):
            d=(X[i] - X[j])*(X[i] - X[j]) + (Y[i] - Y[j])*(Y[i] - Y[j])
            d=np.sqrt(d)/distance
            d_arr=np.append(d_arr,d)
    
           
    DD=np.histogram(d_arr,bins=theta_bins, range=(th_min,th_max))
    N=1.0*n_points*(1.0*n_points-1.0)/2.0
    DD=DD/N
    return DD

def RR_histogram(X,Y,distance,cat_number=1000):
    
    n_points=len(X)
    Xmin=1.0*np.floor( np.amin(X) )
    Xmax=1.0*np.ceil( np.amin(X) )
    Ymin=1.0*np.floor( np.amin(Y) )
    Ymax=1.0*np.ceil( np.amin(Y) )
    

    Xr= Xmin +  ( Xmax - Xmin )*np.random.random_sample(n_points*cat_number)
    Yr=Ymin +   ( Ymax - Ymin )*np.random.random_sample(n_points*cat_number)
    d_arr=[]
    for m in range(cat_number):
        for i in range(n_points):
            for j in range(i,n_points):
                d=(Xr[i + n_points*m ] - Xr[j + n_points*m ])*(Xr[i + n_points*m ] - Xr[j + n_points*m ]) + (Yr[i + n_points*m ] - Yr[j + n_points*m ])*(Yr[i + n_points*m ] - Yr[j + n_points*m ])
                d=np.sqrt(d)/distance
                d_arr=np.append(d_arr,d)
    RR=np.histogram(d_arr,bins=theta_bins, range=(th_min,th_max))
    N=1.0*n_points*(1.0*n_points-1.0)/2.0
    RR=RR/(N*cat_number)            
    return RR

def DR_histogram(X,Y,distance,cat_number=1000):
        n_points=len(X)
    Xmin=1.0*np.floor( np.amin(X) )
    Xmax=1.0*np.ceil( np.amin(X) )
    Ymin=1.0*np.floor( np.amin(Y) )
    Ymax=1.0*np.ceil( np.amin(Y) )
    

    Xr= Xmin +  ( Xmax - Xmin )*np.random.random_sample(n_points*cat_number)
    Yr=Ymin +   ( Ymax - Ymin )*np.random.random_sample(n_points*cat_number)
    d_arr=[]
    for m in range(cat_number): 
            for j in range(i,n_points):
                for i in range(n_points):
                    d=( X[i] - Xr[j + n_points*m] )*( X[i] - Xr[j + n_points*m] ) + ( Y[i] - Yr[j + n_points*m] )*( Y[i] - Yr[j + n_points*m] )
                    d=np.sqrt(d)/distance
                    d_arr=np.append(d_arr,d)
    DR=np.histogram(d_arr,bins=theta_bins, range=(th_min,th_max))
    N=1.0*n_points*(1.0*n_points-1.0)/2.0
    DR=DR/(N*cat_number)            
    return DR


def landy_correlation(DD,RR,DR):
    return (DD - 2.0*DR + RR)/RR
def peebles_correlation(DD,DR):    
    return DD/DR - 1.0
def standard_correlation(DD,RR):
    return DD/RR - 1.0
