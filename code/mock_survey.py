# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

import numpy as np

#The z direction corresponds to the redshift direction
n_x = 5
n_y = 7 
n_z = 6

delta_x = 46.0
delta_y = 35.0


def make_rectangle(width, height):
    vertices = empty((5,2))
    vertices[0,0] = 0.0
    vertices[0,1] = 0.0
    vertices[1,0] = width
    vertices[1,1] = 0.0
    vertices[2,0] = width
    vertices[2,1] = height
    vertices[3,0] = 0.0
    vertices[3,1] = height
    vertices[4,:] = vertices[0,:]
    return vertices


def plot_tile(tiles_x, tiles_y, i_tile, d_x, d_y, color='red'):
    for i_field in range(7):
        rec_base = make_rectangle(delta_x,delta_y)
        plot(rec_base[:,0]+tiles_x[i_tile,i_field]*d_x, rec_base[:,1]+ tiles_y[i_tile,i_field]*d_y, color=color, linewidth=8, alpha=0.7)



#finally dump all the IDs for the fields in the different surveys mimicking the observation
def dump_match(data_path, all_tiles_x, all_tiles_y, all_tiles_z):
    outputfile = data_path+'ID_match_surveys.dat'
    print 'writing to', outputfile
    out = open(outputfile, "w")
    out.write("#ID_survey field_ID i j k\n")
    ID_survey = 0
    for i_z in range(3):    # loop over the slices in Z
        for i_tile in range(5): # loop over the tiles defined in the X-Y plane
        # define a unique ID for each field in the volume
            for i_field in range(7): # these are the fields in the japanese surveys
                index_i = all_tiles_x[i_z, i_tile, i_field]
                index_j = all_tiles_y[i_z, i_tile, i_field]
                index_k = i_z
                field_ID = index_k + n_z*(index_j + (n_y * index_i))
                # print ID_survey, field_ID, 'i,j,k',index_i, index_j, index_k, 'n', n_z, n_y
                out.write("%d %d %d %d %d\n" % (ID_survey, field_ID, index_i, index_j, index_k))
            for i_field in range(3): #tthree correlated fields
                index_i = (all_tiles_x[i_z, i_tile, i_field] + 1)%n_x
                index_j = (all_tiles_y[i_z, i_tile, i_field] + 1 + i_field)%n_y
                index_k = (i_z + 1)%n_z
                field_ID = index_k + n_z*(index_j + (n_y * index_i))
            #print 'corr', ID_survey, field_ID, index_i, index_j, index_k
                out.write("%d %d %d %d %d\n" % (ID_survey, field_ID, index_i, index_j, index_k))
            for i_field in range(2): #two uncorrelated fields
                index_i = (all_tiles_x[i_z, i_tile, i_field] - 1)%n_x
                index_j = (all_tiles_y[i_z, i_tile, i_field] - 1)%n_y
                index_k = (i_z + 2 + i_field)%n_z
                field_ID = index_k + n_z*(index_j + (n_y * index_i))
            #print 'uncorr', ID_survey, field_ID, index_i, index_j, index_k
                out.write("%d %d %d %d %d\n" % (ID_survey, field_ID, index_i, index_j, index_k))
            ID_survey = ID_survey + 1        
    out.close()
    return 
        

# <codecell>

#dump a list of surveys with random sub-boxes but with the same number of fields
def dump_random(data_path):
    outputfile = data_path+'ID_random_surveys.dat'
    print 'writing to', outputfile
    out = open(outputfile, "w")
    out.write("#ID_survey field_ID i j k\n")
    ID_survey = 0
    for i_z in range(3):    # loop over the slices in Z
        for i_tile in range(5): # loop over the tiles defined in the X-Y plane
        # define a unique ID for each field in the volume
            for i_field in range(12): # pick 12 random fields
                index_i = np.random.randint(0,n_x-1)
                index_j = np.random.randint(0,n_y-1)
                index_k = np.random.randint(0,n_z-1)
                field_ID = index_k + n_z*(index_j + (n_y * index_i))
            #print ID_survey, field_ID, index_i, index_j, index_k
                out.write("%d %d %d %d %d\n" % (ID_survey, field_ID, index_i, index_j, index_k))
            ID_survey = ID_survey + 1
    out.close()
    return

#dump 1 single survey with all the sub-boxes
def dump_full(data_path):
    outputfile = data_path+'ID_full_surveys.dat'
    print 'writing to', outputfile
    out = open(outputfile, "w")
    out.write("#ID_survey field_ID i j k\n")
    ID_survey = 0
    for index_k in range(n_z):    # loop over the slices in Z
        for index_j in range(n_y): 
            for  index_i in range(n_x): 
                field_ID = index_k + n_z*(index_j + (n_y * index_i))
            #print ID_survey, field_ID, index_i, index_j, index_k
                out.write("%d %d %d %d %d\n" % (ID_survey, field_ID, index_i, index_j, index_k))
    out.close()
    return 


basic_tiles_x = np.array([[0,1,0,1,0,1,0],[2,3,4,2,3,4,4],[1,0,1,0,1,0,1],[2,3,2,3,4,3,4],[2,2,3,4,2,3,4]])
basic_tiles_y = np.array([[0,0,1,1,2,2,3],[0,0,0,1,1,1,2],[3,4,4,5,5,6,6],[2,2,3,3,3,4,4],[4,5,5,5,6,6,6]])

x_order=np.arange(n_x)
y_order=np.arange(n_y)
z_order=np.arange(n_z)
    
inverted_x = x_order[::-1]
inverted_y = y_order[::-1]

# <codecell>

new_tiles_x_A = np.copy(basic_tiles_x)
new_tiles_y_A = np.copy(basic_tiles_y)
n_tiles = 5
for tile in range(n_tiles):
    new_tiles_x_A[tile,:] = np.copy(basic_tiles_x[tile,:])
    new_tiles_y_A[tile,:] = np.copy(inverted_y[basic_tiles_y[tile,:]])

# <codecell>

new_tiles_x_B = np.copy(basic_tiles_x)
new_tiles_y_B = np.copy(basic_tiles_y)
n_tiles = 5
for tile in range(n_tiles):
    new_tiles_x_B[tile,:] = np.copy(inverted_x[basic_tiles_x[tile,:]])
    new_tiles_y_B[tile,:] = np.copy(basic_tiles_y[tile,:])

# <codecell>

new_tiles_x_C = np.copy(basic_tiles_x)
new_tiles_y_C = np.copy(basic_tiles_y)
n_tiles = 5
for tile in range(n_tiles):
    new_tiles_x_C[tile,:] = np.copy(inverted_x[basic_tiles_x[tile,:]])
    new_tiles_y_C[tile,:] = np.copy(inverted_y[basic_tiles_y[tile,:]])

# <codecell>

new_tiles_x_D = np.copy(basic_tiles_x)
new_tiles_y_D = np.copy(basic_tiles_y)

for i_tile in np.array([0,2]):
    new_tiles_x_D[i_tile,:] = np.copy(basic_tiles_x[i_tile,:] + 3)
    new_tiles_y_D[i_tile,:] = np.copy(basic_tiles_y[i_tile,:])

for i_tile in np.array([1,3,4]):
    new_tiles_x_D[i_tile,:] = np.copy(basic_tiles_x[i_tile,:] - 2)
    new_tiles_y_D[i_tile,:] = np.copy(basic_tiles_y[i_tile,:])
    

# <codecell>

new_tiles_x_E = np.copy(basic_tiles_x)
new_tiles_y_E = np.copy(basic_tiles_y)
n_fields = 5
for field in range(n_fields):
    new_tiles_x_E[field,:] = np.copy(inverted_x[new_tiles_x_D[field,:]])
    new_tiles_y_E[field,:] = np.copy(inverted_y[new_tiles_y_D[field,:]])

# <codecell>

all_tiles_x = np.array([basic_tiles_x, new_tiles_x_A, new_tiles_x_B, new_tiles_x_C, new_tiles_x_D, new_tiles_x_E])
all_tiles_y = np.array([basic_tiles_y, new_tiles_y_A, new_tiles_y_B, new_tiles_y_C, new_tiles_y_D, new_tiles_y_E])
all_tiles_z = np.arange(n_z)

# <codecell>


data_path="../data/mock_survey/"
generate_data = False

if generate_data:
    dump_random(data_path)
    dump_full(data_path)
    dump_match(data_path, all_tiles_x, all_tiles_y, all_tiles_z)



# <codecell>
#i_tile =0
#plot_tile(basic_tiles_x, basic_tiles_y, i_tile, delta_x, delta_y, color='black')
#plot_tile(new_tiles_x_A, new_tiles_y_A, i_tile, delta_x, delta_y, color='blue')
#plot_tile(new_tiles_x_B, new_tiles_y_B, i_tile, delta_x, delta_y, color='green')
#plot_tile(new_tiles_x_C, new_tiles_y_C, i_tile, delta_x, delta_y, color='red')
#plot_tile(new_tiles_x_D, new_tiles_y_D, i_tile, delta_x, delta_y, color='orange')
#plot_tile(new_tiles_x_E, new_tiles_y_E, i_tile, delta_x, delta_y, color='cyan')
