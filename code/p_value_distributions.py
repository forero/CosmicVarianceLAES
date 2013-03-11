# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

from scipy import stats as st
import numpy as np

def ks_test_survey(ID_survey_filename, lae_obs, model_data, filling_fraction=1.0, area=1.0):
    
    #extracts the model data I care about
    field_ID_model = np.int_(model_data[:,0])
    int_i_model = np.int_(model_data[:,1])
    int_j_model = np.int_(model_data[:,2])
    int_k_model = np.int_(model_data[:,3])
    n_halo_model = np.int_(model_data[:,4])
    #field_ID_model = int_k_model + n_z * (int_j_model + (n_y * int_i_model))
  
    #load all data describing the fields that compose the different surveys 
    survey_data = np.loadtxt(ID_survey_filename)
    survey_ID_list = np.int_(survey_data[:,0])
    #print survey_ID_list
    
    #count the total number of surveys included here
    index = np.where(survey_ID_list==survey_ID_list[0])
    n_surveys = survey_data.size/np.size(survey_data[0,:])/np.size(index)

    p_value = np.zeros(n_surveys)
    for survey_ID in range(n_surveys):        
        #select the data for the fields with a given survey ID
        index = np.where(survey_ID_list==survey_ID)
        my_survey = survey_data[index[0],:]
        field_ID_list = np.int_(my_survey[:,1])
        field_i_list = np.int_(my_survey[:,2])
        field_j_list = np.int_(my_survey[:,3])
        field_k_list = np.int_(my_survey[:,4])
    
        #find the halo numbers for each field in all the fields of the survey
        field_nhalo_list = np.zeros(field_ID_list.size)
        for i_field in range(field_ID_list.size):
            index = np.where(field_ID_model==field_ID_list[i_field])
            #print index, field_ID_list[i_field]
            index = index[0]
            #print n_halo_model[index]
            field_nhalo_list[i_field] = n_halo_model[index]
            
        n_lae_theory = filling_fraction * field_nhalo_list/area
        distance,p_value[survey_ID]=st.ks_2samp(n_lae_theory, lae_obs)
    return p_value

# <codecell>
# <codecell>

def get_p_values(full_lae_data_path, full_data_path, survey_filename, filling_fraction, log_min_mass, log_max_mass, mock_area=1.0 ,halo_finder='FOF'):
    filenameout='p_values_'+halo_finder+'_'+survey_filename
    out = open(full_data_path+filenameout, "w")
    out.write("#m_min m_max f_occ survey_ID p_value\n")
    p_values = np.empty((0))
    for i in range(log_min_mass.size):
        for j in range(i,log_max_mass.size):
            m_min = log_min_mass[i]
            m_max = log_max_mass[j]
            filename = 'model_{0}_{1}_{2}.dat'.format(m_min, m_max, 1.0)
            model = np.loadtxt(full_lae_data_path+filename)
            for f_occ in filling_fraction:
                tmp_p_values = ks_test_survey(full_data_path+survey_filename, n_lae_obs, model, filling_fraction=f_occ, area=mock_area)
                n_survey = tmp_p_values.size 
                for survey_ID in range(n_survey):
                    out.write("%f %f %f %d %e\n"%(m_min, m_max, f_occ, survey_ID, tmp_p_values[survey_ID]))
                    p_values = np.append(p_values, tmp_p_values)
    out.close()
         

# <codecell>
# <codecell>



halo_finder = 'FOF'
full_data_path = '../data/mock_survey/'
full_lae_data_path = '../data/laes/'+halo_finder+'/'

filling_fraction = np.arange(0.1,1.1,0.1)
log_min_mass = np.arange(10.0, 12.9, 0.1)
log_max_mass = log_min_mass + 0.1

h=0.7
n_x = 5
n_y = 7
n_z = 6
x_width=46.0
y_width=35.0
D_angular =  1622.4 # in Mpc from Ned Cosmo Calculator
a_exp = 1.0/(3.1 +1.0)
mock_area = (a_exp * x_width/h)*(a_exp * y_width/h) / ((D_angular) * (D_angular))
mock_area = mock_area * (180.0/np.pi) * (180.0/np.pi) * (60.0 * 60.0) #arcmin^2

generate_data = False

if generate_data:

    #load the observational data and convert everything into volumetric densities
    obs_file = "../data/obs/Yamada2012_surface_density.dat"
    obs_prop = np.loadtxt(obs_file, skiprows=4)
    n_lae_obs = obs_prop[:,2]

    survey_filename= "ID_match_surveys.dat"
    get_p_values(full_lae_data_path, full_data_path, survey_filename, filling_fraction, log_min_mass, log_max_mass, mock_area=mock_area, halo_finder=halo_finder)

    survey_filename="ID_random_surveys.dat"
    get_p_values(full_lae_data_path, full_data_path, survey_filename, filling_fraction, log_min_mass, log_max_mass, mock_area=mock_area, halo_finder=halo_finder)

    survey_filename="ID_full_surveys.dat"
    get_p_values(full_lae_data_path, full_data_path, survey_filename, filling_fraction, log_min_mass, log_max_mass, halo_finder=halo_finder)
