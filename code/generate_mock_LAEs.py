#!/usr/bin/env python
# coding: utf-8


import numpy as np
import sys as sys

def generate_all_models():
    
    log_m = np.arange(10.0, 13.1, 0.1)
    box_size = 250.0
    delta_x  = 46.0
    delta_y = 35.0
    delta_z = 41.0
    n_i = int(box_size/delta_x)
    n_j = int(box_size/delta_y)
    n_k = int(box_size/delta_z)
    
    dmh_path="../data/dark_matter/FOF/"
    laes_path="../data/laes/FOF/"

    n_models = n_i * n_j * n_k
    n_bins = log_m.size

    model_numbers = np.zeros((n_i,n_j,n_k,n_bins))

    for i in range(n_i):
        for j in range(n_j):
            for k in range(n_k):
                id0 = k+n_k*(j + n_j*i) 
                print i, j, k
                dmh_filename=dmh_path+"halos_bolshoi_"+str(i)+"-"+str(j)+"-"+str(k)+".csv"
                halos_prop=np.loadtxt(dmh_filename,delimiter=",",skiprows=12) 

                halos_mass=halos_prop[:,4]
                hist, edges = np.histogram(halos_mass, bins=log_m)
                #print hist, edges, hist.size, edges.size, n_bins
                for l in range(n_bins-1):
                    model_numbers[i,j,k,l] = hist[l]
    
    for l in range(n_bins-1):
        for m in range(l+1, n_bins):
            model_name = 'model_{0}_{1}_{2}.dat'.format(log_m[l], log_m[m], 1.0)
            print l, model_name
            out = open(laes_path+model_name, "w")
            for i in range(n_i):
                for j in range(n_j):
                    for k in range(n_k):
                        id0 = k+n_k*(j + n_j*i) 
                        value = sum(model_numbers[i,j,k,l:m])
                        out.write("%d %d %d %d %d\n"%(id0, i, j, k, value))
            out.close()

generate_data = False
if generate_data:
    generate_all_models()



