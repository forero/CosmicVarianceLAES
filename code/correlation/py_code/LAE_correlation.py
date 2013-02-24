import numpy as np 
import correlation_lib as corr


'''--------------------------------------------
survey_type => Indicate the kind of mock survey
to consider for the model.

Description:

match: The mock-field distribution exactly takes the observational distribution
       in position and number. 

random: The number of fields is the same as the observational fields but  
        mock-field positions are taken randomly 
 
full: All the mock-fields are taken into account for the computation

'''
project_path = "/home/jemejia/CosmicVariance/"
def best_model_sel(prob_treshold,survey_type="match", pro_path=project_path):
    ID_file=pro_path + "data/mock_survey/" + "ID_" + survey_type + "_surveys.dat"
    p_values_file= pro_path + "data/mock_survey/" + "p_values_ID_" + survey_type + "_surveys.dat"

    ks_data=np.loadtxt(p_values_file) 
    

    m_min_arr = ks_data[:,0]
    m_max_arr =ks_data[:,1]
    f_occ_arr =ks_data[:,2]
    ID_survey_arr = ks_data[:,3]
    model_prob_arr = ks_data[:,4]
    
    index=np.where(model_prob_array>=prob_treshold)
    best_models=np.empty( [ np.size(index) , np.size( ks_data[0,:] ) ] )
    del(ks_data)

    best_models[:,0]=m_min_arr[index]
    best_models[:,1]=m_max_arr[index]
    best_models[:,2]=f_occ_arr[index]
    best_models[:,3]= ID_survey_arr[index] 
    best_models[:,4]=model_prob_arr[index]

    return best_models
    

best_models=best_model_sel(0.1,survey_type="match", pro_path=project_path )
ID_data=np.loadtxt(ID_file,dtype='int') 
survey_ID=ID_data[:,0]
field_ID=ID_data[:,1]
i_field = ID_data[:,2]
j_field = ID_data[:,3]
j_field = ID_data[:,4]

ID_arr=best_models[:,3]    
for i in range( len(ID_arr) ):
    index=where( survey_ID == ID_arr[i] )
    
